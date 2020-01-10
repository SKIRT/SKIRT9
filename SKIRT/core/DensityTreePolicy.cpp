/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DensityTreePolicy.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

void DensityTreePolicy::setupSelfBefore()
{
    // get the random generator and the number of samples (not a big effort, so we do this even if we don't need it)
    _random = find<Random>();
    _numSamples = find<Configuration>()->numDensitySamples();

    // build lists of media components per material type
    auto ms = find<MediumSystem>(false);  // don't setup the medium system because we are part of it
    int numMedia = ms ? ms->media().size() : 0;
    for (int h = 0; h != numMedia; ++h)
    {
        auto medium = ms->media()[h];
        medium->setup();  // however, we do require each medium component to be setup

        switch (medium->mix()->materialType())
        {
            case MaterialMix::MaterialType::Dust: _dustMedia.push_back(medium); break;
            case MaterialMix::MaterialType::Electrons: _electronMedia.push_back(medium); break;
            case MaterialMix::MaterialType::Gas: _gasMedia.push_back(medium); break;
        }
    }

    // precalculate information for dust
    if (!_dustMedia.empty())
    {
        if (maxDustFraction() > 0) _hasAny = _hasDustAny = _hasDustFraction = true;
        if (maxDustOpticalDepth() > 0) _hasAny = _hasDustAny = _hasDustOpticalDepth = true;
        if (maxDustDensityDispersion() > 0) _hasAny = _hasDustAny = _hasDustDensityDispersion = true;
        if (_hasDustFraction)
        {
            for (auto medium : _dustMedia) _dustMass += medium->mass();
        }
        if (_hasDustOpticalDepth)
        {
            double sigma = 0.;
            double mu = 0.;
            for (auto medium : _dustMedia)
            {
                sigma += medium->mix()->sectionExt(wavelength());
                mu += medium->mix()->mass();
            }
            _dustKappa = sigma / mu;
        }
    }

    // precalculate information for electrons
    if (!_electronMedia.empty() && maxElectronFraction() > 0)
    {
        _hasAny = _hasElectronFraction = true;
        for (auto medium : _electronMedia) _electronNumber += medium->number();
    }

    // precalculate information for gas
    if (!_gasMedia.empty() && maxGasFraction() > 0)
    {
        _hasAny = _hasGasFraction = true;
        for (auto medium : _gasMedia) _gasNumber += medium->number();
    }

    // warn user if none of the criteria were enabled
    if (!_hasAny) find<Log>()->warning("None of the tree subdivision criteria are enabled");
}

////////////////////////////////////////////////////////////////////

bool DensityTreePolicy::needsSubdivide(TreeNode* node)
{
    // results for the sampled mass or number densities, if applicable
    double rho = 0.;          // dust mass density
    double rhomin = DBL_MAX;  // smallest sample for dust mass density
    double rhomax = 0.;       // largest sample for dust mass density
    double ne = 0;            // electron number density
    double ng = 0.;           // gas number density

    // sample densities in node
    if (_hasAny)
    {
        double rhosum = 0;
        double nesum = 0;
        double ngsum = 0;
        for (int i = 0; i != _numSamples; ++i)
        {
            Position bfr = _random->position(node->extent());
            if (_hasDustAny)
            {
                double rhoi = 0.;
                for (auto medium : _dustMedia) rhoi += medium->massDensity(bfr);
                rhosum += rhoi;
                if (rhoi < rhomin) rhomin = rhoi;
                if (rhoi > rhomax) rhomax = rhoi;
            }
            if (_hasElectronFraction)
                for (auto medium : _electronMedia) nesum += medium->numberDensity(bfr);
            if (_hasGasFraction)
                for (auto medium : _gasMedia) ngsum += medium->numberDensity(bfr);
        }
        rho = rhosum / _numSamples;
        ne = nesum / _numSamples;
        ng = ngsum / _numSamples;
    }

    // handle maximum dust mass fraction
    if (_hasDustFraction)
    {
        double delta = rho * node->volume() / _dustMass;
        if (delta > maxDustFraction()) return true;
    }

    // handle maximum dust optical depth
    if (_hasDustOpticalDepth)
    {
        double tau = _dustKappa * rho * node->diagonal();
        if (tau > maxDustOpticalDepth()) return true;
    }

    // handle maximum dust density dispersion
    if (_hasDustDensityDispersion)
    {
        double q = rhomax > 0 ? (rhomax - rhomin) / rhomax : 0.;
        if (q > maxDustDensityDispersion()) return true;
    }

    // handle maximum electron number fraction
    if (_hasElectronFraction)
    {
        double delta = ne * node->volume() / _electronNumber;
        if (delta > maxElectronFraction()) return true;
    }

    // handle maximum gas number fraction
    if (_hasGasFraction)
    {
        double delta = ng * node->volume() / _gasNumber;
        if (delta > maxGasFraction()) return true;
    }

    // if we get here, none of the criteria were violated
    return false;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of nodes evaluated between two invocations of infoIfElapsed()
    const size_t logEvalChunkSize = 10000;

    // maximum number of nodes subdivided between two invocations of infoIfElapsed()
    const size_t logDivideChunkSize = 5000;
}

////////////////////////////////////////////////////////////////////

vector<TreeNode*> DensityTreePolicy::constructTree(TreeNode* root)
{
    auto log = find<Log>();
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // initialize the tree node list with the root node as the first item
    vector<TreeNode*> nodev{root};

    // recursively subdivide the root node until the minimum level has been reached
    int level = 0;    // current level
    size_t lbeg = 0;  // node index range for the current level;
    size_t lend = 1;  // at level 0, the node list contains just the root node
    while (level != minLevel())
    {
        log->info("Subdividing level " + std::to_string(level) + ": " + std::to_string(lend - lbeg) + " nodes");
        log->infoSetElapsed(lend - lbeg);
        for (size_t l = lbeg; l != lend; ++l)
        {
            nodev[l]->subdivide(nodev);
            if ((l + 1) % logDivideChunkSize == 0)
                log->infoIfElapsed("Subdividing level " + std::to_string(level) + ": ", logDivideChunkSize);
        }
        // update iteration variables to the next level
        level++;
        lbeg = lend;
        lend = nodev.size();
    }

    // recursively subdivide the nodes beyond the minimum level until all nodes satisfy the configured criteria
    while (level != maxLevel() && lend != lbeg)
    {
        size_t numEvalNodes = lend - lbeg;
        log->info("Subdividing level " + std::to_string(level) + ": " + std::to_string(numEvalNodes) + " nodes");
        log->infoSetElapsed(numEvalNodes);

        // evaluate nodes at this level: value in the array becomes one for nodes that need to be subdivided
        // we parallelize this operation because it might be resource intensive (e.g. sampling densities)
        Array divide(numEvalNodes);
        parallel->call(numEvalNodes, [this, log, level, lbeg, &nodev, &divide](size_t firstIndex, size_t numIndices) {
            while (numIndices)
            {
                size_t currentChunkSize = min(logEvalChunkSize, numIndices);
                for (size_t l = firstIndex; l != firstIndex + currentChunkSize; ++l)
                {
                    if (needsSubdivide(nodev[lbeg + l])) divide[l] = 1.;
                }
                log->infoIfElapsed("Evaluation for level " + std::to_string(level) + ": ", currentChunkSize);
                firstIndex += currentChunkSize;
                numIndices -= currentChunkSize;
            }
        });
        ProcessManager::sumToAll(divide);

        // subdivide the nodes that have been flagged
        size_t numDivideNodes = divide.sum();
        log->infoSetElapsed(numDivideNodes);
        size_t numDone = 0;
        for (size_t l = 0; l != numEvalNodes; ++l)
        {
            if (divide[l])
            {
                nodev[lbeg + l]->subdivide(nodev);
                numDone++;
                if (numDone % logDivideChunkSize == 0)
                    log->infoIfElapsed("Subdivision for level " + std::to_string(level) + ": ", logDivideChunkSize);
            }
        }
        // update iteration variables to the next level
        level++;
        lbeg = lend;
        lend = nodev.size();
    }

    // sort the neighbors for all nodes
    for (auto node : nodev) node->sortNeighbors();
    return nodev;
}

////////////////////////////////////////////////////////////////////

Range DensityTreePolicy::wavelengthRange() const
{
    if (maxDustOpticalDepth() > 0)
        return Range(wavelength(), wavelength());
    else
        return Range();
}

////////////////////////////////////////////////////////////////////

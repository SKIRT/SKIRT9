/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DensityTreePolicy.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MassInBoxInterface.hpp"
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
    TreePolicy::setupSelfBefore();

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

    // build corresponding lists of MassInBoxInterface pointers
    auto buildMIBlist = [](vector<MassInBoxInterface*>& mibv, const vector<Medium*>& media) {
        for (auto medium : media)
        {
            auto mib = medium->interface<MassInBoxInterface>(0, false);
            if (!mib)
            {
                mibv.clear();  // if any medium component lacks support, clear the complete list
                return;
            }
            mibv.push_back(mib);
        }
    };
    buildMIBlist(_dustMIBv, _dustMedia);
    buildMIBlist(_electronMIBv, _electronMedia);
    buildMIBlist(_gasMIBv, _gasMedia);

    // precalculate information for dust
    if (!_dustMedia.empty())
    {
        if (maxDustFraction() > 0) _hasDustFraction = true;
        if (maxDustOpticalDepth() > 0) _hasDustOpticalDepth = true;
        if (maxDustDensityDispersion() > 0) _hasDustDensityDispersion = true;
        _needDustSamples = _dustMIBv.empty() | _hasDustDensityDispersion;
        _hasDustMIB = !_dustMIBv.empty();
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
        _hasElectronFraction = true;
        _needElectronSamples = _electronMIBv.empty();
        _hasElectronMIB = !_electronMIBv.empty();
        for (auto medium : _electronMedia) _electronNumber += medium->number();
    }

    // precalculate information for gas
    if (!_gasMedia.empty() && maxGasFraction() > 0)
    {
        _hasGasFraction = true;
        _needGasSamples = _gasMIBv.empty();
        _hasGasMIB = !_gasMIBv.empty();
        for (auto medium : _gasMedia) _gasNumber += medium->number();
    }

    // determine composite flags
    _hasAny =
        _hasDustFraction | _hasDustOpticalDepth | _hasDustDensityDispersion | _hasElectronFraction | _hasGasFraction;
    _needAnySamples = _needDustSamples | _needElectronSamples | _needGasSamples;

    // warn user if none of the criteria were enabled
    if (!_hasAny) find<Log>()->warning("None of the tree subdivision criteria are enabled");
}

////////////////////////////////////////////////////////////////////

bool DensityTreePolicy::needsSubdivide(TreeNode* node)
{
    // handle minimum and maximum level
    if (node->level() < minLevel()) return true;
    if (node->level() >= maxLevel()) return false;

    // results for the sampled mass or number densities, if applicable
    double rho = 0.;          // dust mass density
    double rhomin = DBL_MAX;  // smallest sample for dust mass density
    double rhomax = 0.;       // largest sample for dust mass density
    double ne = 0;            // electron number density
    double ng = 0.;           // gas number density

    // sample densities in node, if needed
    if (_needAnySamples)
    {
        double rhosum = 0;
        double nesum = 0;
        double ngsum = 0;
        for (int i = 0; i != _numSamples; ++i)
        {
            Position bfr = _random->position(node->extent());
            if (_needDustSamples)
            {
                double rhoi = 0.;
                for (auto medium : _dustMedia) rhoi += medium->massDensity(bfr);
                rhosum += rhoi;
                if (rhoi < rhomin) rhomin = rhoi;
                if (rhoi > rhomax) rhomax = rhoi;
            }
            if (_needElectronSamples)
                for (auto medium : _electronMedia) nesum += medium->numberDensity(bfr);
            if (_needGasSamples)
                for (auto medium : _gasMedia) ngsum += medium->numberDensity(bfr);
        }
        rho = rhosum / _numSamples;
        ne = nesum / _numSamples;
        ng = ngsum / _numSamples;
    }

    // results for the volume-integrated mass or number
    double M = 0.;   // volume-integrated dust mass
    double Ne = 0.;  // volume-integrated electron number
    double Ng = 0.;  // volume-integrated gas number

    double V = node->volume();  // node volume

    // calculate volume-integrated quantities and densities from MIB if possible, and else from samples
    if (_hasDustMIB)
    {
        for (auto mib : _dustMIBv) M += mib->massInBox(node->extent());
        rho = M / V;
    }
    else
    {
        M = rho * V;
    }
    if (_hasElectronMIB)
    {
        for (auto mib : _electronMIBv) Ne += mib->numberInBox(node->extent());
    }
    else
    {
        Ne = ne * V;
    }
    if (_hasGasMIB)
    {
        for (auto mib : _gasMIBv) Ng += mib->numberInBox(node->extent());
    }
    else
    {
        Ng = ng * V;
    }

    // handle maximum dust mass fraction
    if (_hasDustFraction)
    {
        double delta = M / _dustMass;
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
        double delta = Ne / _electronNumber;
        if (delta > maxElectronFraction()) return true;
    }

    // handle maximum gas number fraction
    if (_hasGasFraction)
    {
        double delta = Ng / _gasNumber;
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

    // inform user about the use of MassInBoxInterface
    if (_hasDustMIB) find<Log>()->info("  (obtaining dust densities through calculation rather than sampling)");
    if (_hasElectronMIB) find<Log>()->info("  (obtaining electron densities through calculation rather than sampling)");
    if (_hasGasMIB) find<Log>()->info("  (obtaining gas densities through calculation rather than sampling)");

    // initialize the tree node list with the root node as the first item
    vector<TreeNode*> nodev{root};

    // initialize iteration variables to level 0
    int level = 0;    // current level
    size_t lbeg = 0;  // node index range for the current level;
    size_t lend = 1;  // at level 0, the node list contains just the root node

    // recursively subdivide nodes until all nodes satisfy the configured criteria
    while (lend != lbeg)
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

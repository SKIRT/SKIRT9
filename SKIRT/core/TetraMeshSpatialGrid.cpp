/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TetraMeshSpatialGrid.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SiteListInterface.hpp"
#include "TetraMeshInterface.hpp"
#include "TetraMeshSnapshot.hpp"

//////////////////////////////////////////////////////////////////////

TetraMeshSpatialGrid::~TetraMeshSpatialGrid() {}

//////////////////////////////////////////////////////////////////////

namespace
{
    // sample random positions from the given media components with the given relative weigths
    vector<Vec> sampleMedia(const vector<Medium*>& media, const vector<double>& weights, const Box& extent,
                            int numSites)
    {
        if (media.empty()) throw FATALERROR("There is no material of the requested type in the medium system");

        // build the cumulative weight distribution
        Array Xv;
        NR::cdf(Xv, media.size(), [weights](int h) { return weights[h]; });
        auto random = media[0]->find<Random>();

        // sample
        vector<Vec> rv(numSites);
        for (int m = 0; m != numSites;)
        {
            int h = NR::locateClip(Xv, random->uniform());
            Position p = media[h]->generatePosition();
            if (extent.contains(p)) rv[m++] = p;  // discard any points outside of the domain
        }
        return rv;
    }
}

//////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::setupSelfBefore()
{
    BoxSpatialGrid::setupSelfBefore();

    _random = find<Random>();
    _numSamples = find<Configuration>()->numDensitySamples();

    const MediumSystem* ms = find<MediumSystem>();
    if (!ms) throw FATALERROR("MediumSystem not found");

    for (Medium* medium : ms->media())
    {
        medium->setup();
        switch (medium->mix()->materialType())
        {
            case MaterialMix::MaterialType::Dust: _dustMedia.push_back(medium); break;
            case MaterialMix::MaterialType::Electrons: _electronMedia.push_back(medium); break;
            case MaterialMix::MaterialType::Gas: _gasMedia.push_back(medium); break;
        }
    }

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

    _mesh = new TetraMeshSnapshot(this, extent());
}

bool TetraMeshSpatialGrid::tetUnsuitable(double* pa, double* pb, double* pc, double* pd, double vol) const
{
    Vec a(pa[0], pa[1], pa[2]);
    Vec b(pb[0], pb[1], pb[2]);
    Vec c(pc[0], pc[1], pc[2]);
    Vec d(pd[0], pd[1], pd[2]);

    // results for the sampled mass or number densities, if applicable
    double rho = 0.;  // dust mass density
    // double rhovar = 0.;       // dust mass variance
    double rhomin = DBL_MAX;  // smallest sample for dust mass density
    double rhomax = 0.;       // largest sample for dust mass density
    double ne = 0;            // electron number density
    double ng = 0.;           // gas number density

    // sample densities in node
    if (_hasAny)
    {
        double rhosum = 0;
        // double rhosum2 = 0;
        double nesum = 0;
        double ngsum = 0;
        for (int i = 0; i != _numSamples; ++i)
        {
            double s = random()->uniform();
            double t = random()->uniform();
            double u = random()->uniform();
            double r = Tetra::generateBarycentric(s, t, u);
            Position bfr(r * a + u * b + t * c + s * d);
            if (_hasDustAny)
            {
                double rhoi = 0.;
                for (auto medium : _dustMedia) rhoi += medium->massDensity(bfr);
                rhosum += rhoi;
                // rhosum2 += rhoi * rhoi;
                if (rhoi < rhomin) rhomin = rhoi;
                if (rhoi > rhomax) rhomax = rhoi;
            }
            if (_hasElectronFraction)
                for (auto medium : _electronMedia) nesum += medium->numberDensity(bfr);
            if (_hasGasFraction)
                for (auto medium : _gasMedia) ngsum += medium->numberDensity(bfr);
        }
        rho = rhosum / _numSamples;
        // rhovar = rhomax > 0 ? (rhosum2 / _numSamples - rho * rho) * (vol / _dustMass) * (vol / _dustMass) / (maxDustFraction() * maxDustFraction()) : 0.;
        ne = nesum / _numSamples;
        ng = ngsum / _numSamples;
    }

    // handle maximum dust mass fraction
    if (_hasDustFraction)
    {
        double delta = rho * vol / _dustMass;
        if (delta > maxDustFraction()) return true;
        if (delta < minDustFraction()) return false;
    }

    // handle maximum dust optical depth
    // if (_hasDustOpticalDepth)
    // {
    //     double tau = _dustKappa * rho * node->diagonal();
    //     if (tau > maxDustOpticalDepth()) return true;
    // }

    // handle maximum dust density dispersion
    if (_hasDustDensityDispersion)
    {
        // if (rhovar > maxDustDensityDispersion()) return true;

        double q = rhomax > 0 ? (rhomax - rhomin) / rhomax : 0.;
        if (q > maxDustDensityDispersion()) return true;
    }

    // handle maximum electron number fraction
    if (_hasElectronFraction)
    {
        double delta = ne * vol / _electronNumber;
        if (delta > maxElectronFraction()) return true;
    }

    // handle maximum gas number fraction
    if (_hasGasFraction)
    {
        double delta = ng * vol / _gasNumber;
        if (delta > maxGasFraction()) return true;
    }

    // if we get here, none of the criteria were violated
    return false;
}

//////////////////////////////////////////////////////////////////////

int TetraMeshSpatialGrid::numCells() const
{
    return _mesh->numEntities();
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::volume(int m) const
{
    return _mesh->volume(m);
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::diagonal(int m) const
{
    return cbrt(3. * _mesh->volume(m));
}

//////////////////////////////////////////////////////////////////////

int TetraMeshSpatialGrid::cellIndex(Position bfr) const
{
    return _mesh->cellIndex(bfr);
}

//////////////////////////////////////////////////////////////////////

Position TetraMeshSpatialGrid::centralPositionInCell(int m) const
{
    return _mesh->position(m);
}

//////////////////////////////////////////////////////////////////////

Position TetraMeshSpatialGrid::randomPositionInCell(int m) const
{
    return _mesh->generatePosition(m);
}

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> TetraMeshSpatialGrid::createPathSegmentGenerator() const
{
    return _mesh->createPathSegmentGenerator();
}

//////////////////////////////////////////////////////////////////////

void TetraMeshSpatialGrid::writeGridPlotFiles(const SimulationItem* probe) const
{
    _mesh->writeGridPlotFiles(probe);
}

//////////////////////////////////////////////////////////////////////

double TetraMeshSpatialGrid::numberDensity(int /*h*/, int m) const
{
    return _norm * _mesh->density(m);
}

//////////////////////////////////////////////////////////////////////

bool TetraMeshSpatialGrid::offersInterface(const std::type_info& interfaceTypeInfo) const
{
    if (interfaceTypeInfo == typeid(DensityInCellInterface))
    {
        return false;
        // if (_policy != Policy::ImportedMesh) return false;
        // auto ms = find<MediumSystem>(false);
        // return ms && ms->media().size() == 1;
    }
    return SpatialGrid::offersInterface(interfaceTypeInfo);
}

//////////////////////////////////////////////////////////////////////

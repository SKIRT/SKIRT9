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

    // determine an appropriate set of sites and construct the Tetra mesh
    switch (_policy)
    {
        case Policy::Uniform:
        {
            auto random = find<Random>();
            vector<Vec> rv(_numSites);
            for (int m = 0; m != _numSites; ++m) rv[m] = random->position(extent());
            _mesh = new TetraMeshSnapshot(this, extent(), rv, _mindihedral);
            break;
        }
        case Policy::CentralPeak:
        {
            auto random = find<Random>();
            const int a = 1000;  // steepness of the peak; the central 1/a portion is NOT covered
            const double rscale = extent().rmax().norm();
            vector<Vec> rv(_numSites);
            for (int m = 1; m != _numSites;)  // skip first particle so that it remains (0,0,0)
            {
                double r = rscale * pow(1. / a, random->uniform());  // random distribution according to 1/x
                Direction k = random->direction();
                Position p = Position(r, k);
                if (extent().contains(p)) rv[m++] = p;  // discard any points outside of the domain
            }
            _mesh = new TetraMeshSnapshot(this, extent(), rv, _mindihedral);
            break;
        }
        case Policy::DustDensity:
        {
            // build a list of media that have this material type with corresponding weights
            vector<Medium*> media;
            vector<double> weights;
            auto ms = find<MediumSystem>();
            for (auto medium : ms->media())
                if (medium->mix()->isDust()) media.push_back(medium);
            for (auto medium : media) weights.push_back(medium->mass());
            _mesh =
                new TetraMeshSnapshot(this, extent(), sampleMedia(media, weights, extent(), _numSites), _mindihedral);
            break;
        }
        case Policy::ElectronDensity:
        {
            // build a list of media that have this material type with corresponding weights
            vector<Medium*> media;
            vector<double> weights;
            auto ms = find<MediumSystem>();
            for (auto medium : ms->media())
                if (medium->mix()->isElectrons()) media.push_back(medium);
            for (auto medium : media) weights.push_back(medium->number());
            _mesh =
                new TetraMeshSnapshot(this, extent(), sampleMedia(media, weights, extent(), _numSites), _mindihedral);
            break;
        }
        case Policy::GasDensity:
        {
            // build a list of media that have this material type with corresponding weights
            vector<Medium*> media;
            vector<double> weights;
            auto ms = find<MediumSystem>();
            for (auto medium : ms->media())
                if (medium->mix()->isGas()) media.push_back(medium);
            for (auto medium : media) weights.push_back(medium->number());
            _mesh =
                new TetraMeshSnapshot(this, extent(), sampleMedia(media, weights, extent(), _numSites), _mindihedral);
            break;
        }
        case Policy::File:
        {
            _mesh = new TetraMeshSnapshot(this, extent(), _filename, _mindihedral);
            break;
        }
        case Policy::ImportedSites:
        {
            auto sli = find<MediumSystem>()->interface<SiteListInterface>(2);
            _mesh = new TetraMeshSnapshot(this, extent(), sli, _mindihedral);
            break;
        }
        case Policy::ImportedMesh:
        {
            auto ms = find<MediumSystem>(false);
            _mesh = ms->interface<TetraMeshInterface>(2)->tetraMesh();

            // if there is a single medium component, calculate the normalization factor imposed by it;
            // we need this to directly compute cell densities for the DensityInCellInterface
            if (ms->media().size() == 1) _norm = _mesh->mass() > 0 ? ms->media()[0]->number() / _mesh->mass() : 0.;
            break;
        }
    }
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

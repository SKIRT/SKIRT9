/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CellSnapshot.hpp"
#include "EntityCollection.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

Box CellSnapshot::boxForCell(int m) const
{
    const auto& prop = _propv[m];
    int i = boxIndex();
    return Box(prop[i], prop[i + 1], prop[i + 2], prop[i + 3], prop[i + 4], prop[i + 5]);
}

////////////////////////////////////////////////////////////////////

void CellSnapshot::readAndClose()
{
    // read the snapshot cell info into memory
    _propv = infile()->readAllRows();

    // close the file
    close();

    // inform the user
    log()->info("  Number of cells: " + std::to_string(_propv.size()));

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy()) calculateDensityAndMass(_rhov, _cumrhov, _mass);

    // if needed, construct a search structure for the cells
    if (hasMassDensityPolicy() || needGetEntities())
    {
        log()->info("Constructing search grid for " + std::to_string(_propv.size()) + " cells...");
        auto bounds = [this](int m) { return boxForCell(m); };
        auto intersects = [this](int m, const Box& box) { return box.intersects(boxForCell(m)); };
        _search.loadEntities(_propv.size(), bounds, intersects);

        int nb = _search.numBlocks();
        log()->info("  Number of blocks in grid: " + std::to_string(nb * nb * nb) + " (" + std::to_string(nb) + "^3)");
        log()->info("  Smallest number of cells per block: " + std::to_string(_search.minEntitiesPerBlock()));
        log()->info("  Largest  number of cells per block: " + std::to_string(_search.maxEntitiesPerBlock()));
        log()->info("  Average  number of cells per block: "
                    + StringUtils::toString(_search.avgEntitiesPerBlock(), 'f', 1));
    }
}

////////////////////////////////////////////////////////////////////

Box CellSnapshot::extent() const
{
    // if there are no cells, return an empty box
    if (_propv.empty()) return Box();

    // if there is a search structure, ask it to return the extent (it is already calculated)
    if (_search.numBlocks()) return _search.extent();

    // otherwise find the spatial range of the cells
    double xmin = +std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = +std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmin = +std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (const Array& prop : _propv)
    {
        xmin = min(xmin, prop[boxIndex() + 0]);
        xmax = max(xmax, prop[boxIndex() + 3]);
        ymin = min(ymin, prop[boxIndex() + 1]);
        ymax = max(ymax, prop[boxIndex() + 4]);
        zmin = min(zmin, prop[boxIndex() + 2]);
        zmax = max(zmax, prop[boxIndex() + 5]);
    }
    return Box(xmin, ymin, zmin, xmax, ymax, zmax);
}

////////////////////////////////////////////////////////////////////

int CellSnapshot::numEntities() const
{
    return _propv.size();
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::volume(int m) const
{
    return boxForCell(m).volume();
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::density(Position bfr) const
{
    for (int m : _search.entitiesFor(bfr))
    {
        if (boxForCell(m).contains(bfr)) return _rhov[m];
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double CellSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position CellSnapshot::position(int m) const
{
    return Position(boxForCell(m).center());
}

////////////////////////////////////////////////////////////////////

Position CellSnapshot::generatePosition(int m) const
{
    return random()->position(boxForCell(m));
}

////////////////////////////////////////////////////////////////////

Position CellSnapshot::generatePosition() const
{
    // if there are no cells, return the origin
    if (_propv.empty()) return Position();

    // select a cell according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

const Array& CellSnapshot::properties(int m) const
{
    return _propv[m];
}

////////////////////////////////////////////////////////////////////

void CellSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    for (int m : _search.entitiesFor(bfr))
    {
        if (boxForCell(m).contains(bfr))
        {
            entities.addSingle(m);
            return;
        }
    }
    entities.clear();
}

////////////////////////////////////////////////////////////////////

void CellSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
{
    entities.clear();
    for (int m : _search.entitiesFor(bfr, bfk))
    {
        double smin, smax;
        if (boxForCell(m).intersects(bfr, bfk, smin, smax)) entities.add(m, smax - smin);
    };
}

////////////////////////////////////////////////////////////////////

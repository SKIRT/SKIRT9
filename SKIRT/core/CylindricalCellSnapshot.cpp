/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylindricalCellSnapshot.hpp"
#include "EntityCollection.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void CylindricalCellSnapshot::readAndClose()
{
    // read the snapshot cell info into memory
    _propv = infile()->readAllRows();

    // close the file
    close();

    // inform the user
    log()->info("  Number of cells: " + std::to_string(_propv.size()));

    // build CylCell objects so we can calculate volume and other geometric properties
    _cellv.reserve(_propv.size());
    int i = boxIndex();
    for (const Array& prop : _propv)
        _cellv.emplace_back(prop[i], prop[i + 1], prop[i + 2], prop[i + 3], prop[i + 4], prop[i + 5]);

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy()) calculateDensityAndMass(_rhov, _cumrhov, _mass);

    // if needed, construct a search structure for the cells
    if (hasMassDensityPolicy() || needGetEntities())
    {
        log()->info("Constructing search grid for " + std::to_string(_propv.size()) + " cells...");
        auto bounds = [this](int m) { return _cellv[m].boundingBox(); };
        auto intersects = [this](int m, const Box& box) { return box.intersects(_cellv[m].boundingBox()); };
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

Box CylindricalCellSnapshot::extent() const
{
    // if there is a search structure, ask it to return the extent (it is already calculated)
    if (_search.numBlocks()) return _search.extent();

    // otherwise find the spatial range of the cells
    Box extent;
    for (const auto& cell : _cellv) extent.extend(cell.boundingBox());
    return extent;
}

////////////////////////////////////////////////////////////////////

int CylindricalCellSnapshot::numEntities() const
{
    return _propv.size();
}

////////////////////////////////////////////////////////////////////

double CylindricalCellSnapshot::volume(int m) const
{
    return _cellv[m].volume();
}

////////////////////////////////////////////////////////////////////

double CylindricalCellSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double CylindricalCellSnapshot::density(Position bfr) const
{
    for (int m : _search.entitiesFor(bfr))
    {
        if (_cellv[m].contains(bfr)) return _rhov[m];
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double CylindricalCellSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position CylindricalCellSnapshot::position(int m) const
{
    return Position(_cellv[m].center());
}

////////////////////////////////////////////////////////////////////

Position CylindricalCellSnapshot::generatePosition(int m) const
{
    double Rmin, phimin, zmin, Rmax, phimax, zmax;
    _cellv[m].extent(Rmin, phimin, zmin, Rmax, phimax, zmax);

    double R = sqrt(Rmin * Rmin + (Rmax - Rmin) * (Rmax + Rmin) * random()->uniform());
    double phi = phimin + (phimax - phimin) * random()->uniform();
    double z = zmin + (zmax - zmin) * random()->uniform();
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

////////////////////////////////////////////////////////////////////

Position CylindricalCellSnapshot::generatePosition() const
{
    // if there are no cells, return the origin
    if (_propv.empty()) return Position();

    // select a cell according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

const Array& CylindricalCellSnapshot::properties(int m) const
{
    return _propv[m];
}

////////////////////////////////////////////////////////////////////

void CylindricalCellSnapshot::getEntities(EntityCollection& entities, Position bfr) const
{
    for (int m : _search.entitiesFor(bfr))
    {
        if (_cellv[m].contains(bfr))
        {
            entities.addSingle(m);
            return;
        }
    }
    entities.clear();
}

////////////////////////////////////////////////////////////////////

void CylindricalCellSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
{
    entities.clear();
    for (int m : _search.entitiesFor(bfr, bfk))
    {
        entities.add(m, _cellv[m].intersection(bfr, bfk));
    };
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CylindricalCellSnapshot.hpp"
#include "EntityCollection.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void CylindricalCellSnapshot::setNumAutoRevolveBins(int numBins)
{
    _numAutoRevolveBins = numBins;
}

////////////////////////////////////////////////////////////////////

void CylindricalCellSnapshot::readAndClose()
{
    // read the snapshot cell info into memory
    _propv = infile()->readAllRows();

    // close the file
    close();

    // perform auto-revolve if requested
    if (_numAutoRevolveBins >= 2)
    {
        int num2DCells = _propv.size();
        log()->info("  Auto-revolving " + std::to_string(num2DCells) + " 2D cells using "
                    + std::to_string(_numAutoRevolveBins) + " azimuth bins");

        // construct the azimuth grid, and make sure the endpoints are exact
        Array phiv;
        NR::buildLinearGrid(phiv, -M_PI, M_PI, _numAutoRevolveBins);
        phiv[0] = -M_PI;
        phiv[_numAutoRevolveBins] = M_PI;

        // get indices for columns we need during the revolve operation (mass indices might be -1)
        int iphimin = boxIndex() + 1;
        int iphimax = boxIndex() + 4;
        int imass1 = massIndex();
        int imass2 = initialMassIndex();
        int imass3 = currentMassIndex();

        // move the original 2D cells to a temporary vector
        vector<Array> prop2Dv;
        _propv.swap(prop2Dv);
        _propv.reserve(num2DCells * _numAutoRevolveBins);

        // loop over the original 2D cells
        for (Array& prop : prop2Dv)
        {
            // verify that the cell is 2D
            if (prop[iphimin] || prop[iphimax]) throw FATALERROR("2D cell in input file has nonzero azimuth angle(s)");

            // distribute masses over azimuth bins
            if (imass1 >= 0) prop[imass1] /= _numAutoRevolveBins;
            if (imass2 >= 0) prop[imass2] /= _numAutoRevolveBins;
            if (imass3 >= 0) prop[imass3] /= _numAutoRevolveBins;

            // loop over azimuth bins and add a new 3D cell for each
            for (int k = 0; k != _numAutoRevolveBins; ++k)
            {
                prop[iphimin] = phiv[k];
                prop[iphimax] = phiv[k + 1];
                _propv.emplace_back(prop);
            }
        }
    }

    // inform the user
    log()->info("  Number of cells: " + std::to_string(_propv.size()));

    // build CylindricalCell objects so we can calculate volume and other geometric properties
    _cellv.reserve(_propv.size());
    int bi = boxIndex();
    for (const Array& prop : _propv)
        _cellv.emplace_back(prop[bi], prop[bi + 1], prop[bi + 2], prop[bi + 3], prop[bi + 4], prop[bi + 5]);

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy()) calculateDensityAndMass(_rhov, _cumrhov, _mass);

    // if applicable, convert velocity and/or magnetic field vectors from cylindrical to Cartesian coordinates
    int vi = velocityIndex();
    int mi = magneticFieldIndex();
    if (vi >= 0 || mi >= 0)
    {
        for (Array& prop : _propv)
        {
            // get the central angle of the cell
            double phi = 0.5 * (prop[bi + 1] + prop[bi + 4]);
            double sinphi = sin(phi);
            double cosphi = cos(phi);

            // convert the velocity vector
            if (vi >= 0)
            {
                double vR = prop[vi];
                double vphi = prop[vi + 1];
                prop[vi] = vR * cosphi - vphi * sinphi;
                prop[vi + 1] = vR * sinphi + vphi * cosphi;
            }

            // convert the magnetic field vector
            if (mi >= 0)
            {
                double BR = prop[mi];
                double Bphi = prop[mi + 1];
                prop[mi] = BR * cosphi - Bphi * sinphi;
                prop[mi + 1] = BR * sinphi + Bphi * cosphi;
            }
        }
    }

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
    return _cellv[m].center();
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

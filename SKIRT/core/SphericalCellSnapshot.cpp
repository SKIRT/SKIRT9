/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalCellSnapshot.hpp"
#include "Cubic.hpp"
#include "EntityCollection.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void SphericalCellSnapshot::setNumAutoRevolveBins(int numInclinationBins, int numAzimuthBins)
{
    _numAutoInclinationBins = numInclinationBins;
    _numAutoAzimuthBins = numAzimuthBins;
}

////////////////////////////////////////////////////////////////////

void SphericalCellSnapshot::readAndClose()
{
    // read the snapshot cell info into memory
    _propv = infile()->readAllRows();

    // close the file
    close();

    // perform inclination auto-revolve if requested
    if (_numAutoInclinationBins >= 2)
    {
        int numCells = _propv.size();
        log()->info("  Auto-revolving " + std::to_string(numCells) + " cells using "
                    + std::to_string(_numAutoInclinationBins) + " inclination bins");

        // construct the inclination grid, and make sure the endpoints are exact;
        // use cosine grid point placement so that all inclination bins have the same volume
        Array thetav;
        NR::buildLinearGrid(thetav, 1., -1., _numAutoInclinationBins);
        thetav = acos(thetav);
        thetav[0] = 0.;
        thetav[_numAutoInclinationBins] = M_PI;

        // get indices for columns we need during the revolve operation (mass indices might be -1)
        int ithetamin = boxIndex() + 1;
        int ithetamax = boxIndex() + 4;
        int imass1 = massIndex();
        int imass2 = initialMassIndex();
        int imass3 = currentMassIndex();

        // move the original cells to a temporary vector
        vector<Array> propOrigv;
        _propv.swap(propOrigv);
        _propv.reserve(numCells * _numAutoInclinationBins);

        // loop over the original cells
        for (Array& prop : propOrigv)
        {
            // verify that the cell is non-3D
            if (prop[ithetamin] || prop[ithetamax])
                throw FATALERROR("non-3D cell in input file has nonzero inclination angle(s)");

            // distribute masses over inclination bins
            if (imass1 >= 0) prop[imass1] /= _numAutoInclinationBins;
            if (imass2 >= 0) prop[imass2] /= _numAutoInclinationBins;
            if (imass3 >= 0) prop[imass3] /= _numAutoInclinationBins;

            // loop over inclination bins and add a new 3D cell for each
            for (int j = 0; j != _numAutoInclinationBins; ++j)
            {
                prop[ithetamin] = thetav[j];
                prop[ithetamax] = thetav[j + 1];
                _propv.push_back(prop);
            }
        }
    }

    // perform azimuth auto-revolve if requested
    if (_numAutoAzimuthBins >= 2)
    {
        int numCells = _propv.size();
        log()->info("  Auto-revolving " + std::to_string(numCells) + " cells using "
                    + std::to_string(_numAutoAzimuthBins) + " azimuth bins");

        // construct the azimuth grid, and make sure the endpoints are exact
        Array phiv;
        NR::buildLinearGrid(phiv, -M_PI, M_PI, _numAutoAzimuthBins);
        phiv[0] = -M_PI;
        phiv[_numAutoAzimuthBins] = M_PI;

        // get indices for columns we need during the revolve operation (mass indices might be -1)
        int iphimin = boxIndex() + 2;
        int iphimax = boxIndex() + 5;
        int imass1 = massIndex();
        int imass2 = initialMassIndex();
        int imass3 = currentMassIndex();

        // move the original cells to a temporary vector
        vector<Array> propOrigv;
        _propv.swap(propOrigv);
        _propv.reserve(numCells * _numAutoAzimuthBins);

        // loop over the original cells
        for (Array& prop : propOrigv)
        {
            // verify that the cell is non-3D
            if (prop[iphimin] || prop[iphimax])
                throw FATALERROR("non-3D cell in input file has nonzero azimuth angle(s)");

            // distribute masses over azimuth bins
            if (imass1 >= 0) prop[imass1] /= _numAutoAzimuthBins;
            if (imass2 >= 0) prop[imass2] /= _numAutoAzimuthBins;
            if (imass3 >= 0) prop[imass3] /= _numAutoAzimuthBins;

            // loop over azimuth bins and add a new 3D cell for each
            for (int k = 0; k != _numAutoAzimuthBins; ++k)
            {
                prop[iphimin] = phiv[k];
                prop[iphimax] = phiv[k + 1];
                _propv.push_back(prop);
            }
        }
    }

    // inform the user
    log()->info("  Number of cells: " + std::to_string(_propv.size()));

    // build SphericalCell objects so we can calculate volume and other geometric properties
    _cellv.reserve(_propv.size());
    int bi = boxIndex();
    for (const Array& prop : _propv)
        _cellv.emplace_back(prop[bi], prop[bi + 1], prop[bi + 2], prop[bi + 3], prop[bi + 4], prop[bi + 5]);

    // if a mass density policy has been set, calculate masses and densities for all cells
    if (hasMassDensityPolicy()) calculateDensityAndMass(_rhov, _cumrhov, _mass);

    // if applicable, convert velocity and/or magnetic field vectors from spherical to Cartesian coordinates
    int vi = velocityIndex();
    int mi = magneticFieldIndex();
    if (vi >= 0 || mi >= 0)
    {
        for (Array& prop : _propv)
        {
            // get the central inclination angle of the cell
            double theta = 0.5 * (prop[bi + 1] + prop[bi + 4]);
            double sintheta = sin(theta);
            double costheta = cos(theta);

            // get the central azimuth angle of the cell
            double phi = 0.5 * (prop[bi + 2] + prop[bi + 5]);
            double sinphi = sin(phi);
            double cosphi = cos(phi);

            // convert the velocity vector
            if (vi >= 0)
            {
                double vr = prop[vi];
                double vtheta = prop[vi + 1];
                double vphi = prop[vi + 2];
                prop[vi + 0] = vr * sintheta * cosphi + vtheta * costheta * cosphi - vphi * sinphi;
                prop[vi + 1] = vr * sintheta * sinphi + vtheta * costheta * sinphi + vphi * cosphi;
                prop[vi + 2] = vr * costheta - vtheta * sintheta;
            }

            // convert the magnetic field vector
            if (mi >= 0)
            {
                double Br = prop[mi];
                double Btheta = prop[mi + 1];
                double Bphi = prop[mi + 2];
                prop[mi + 0] = Br * sintheta * cosphi + Btheta * costheta * cosphi - Bphi * sinphi;
                prop[mi + 1] = Br * sintheta * sinphi + Btheta * costheta * sinphi + Bphi * cosphi;
                prop[mi + 2] = Br * costheta - Btheta * sintheta;
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

Box SphericalCellSnapshot::extent() const
{
    // if there is a search structure, ask it to return the extent (it is already calculated)
    if (_search.numBlocks()) return _search.extent();

    // otherwise find the spatial range of the cells
    Box extent;
    for (const auto& cell : _cellv) extent.extend(cell.boundingBox());
    return extent;
}

////////////////////////////////////////////////////////////////////

int SphericalCellSnapshot::numEntities() const
{
    return _propv.size();
}

////////////////////////////////////////////////////////////////////

double SphericalCellSnapshot::volume(int m) const
{
    return _cellv[m].volume();
}

////////////////////////////////////////////////////////////////////

double SphericalCellSnapshot::density(int m) const
{
    return _rhov[m];
}

////////////////////////////////////////////////////////////////////

double SphericalCellSnapshot::density(Position bfr) const
{
    for (int m : _search.entitiesFor(bfr))
    {
        if (_cellv[m].contains(bfr)) return _rhov[m];
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SphericalCellSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position SphericalCellSnapshot::position(int m) const
{
    return _cellv[m].center();
}

////////////////////////////////////////////////////////////////////

Position SphericalCellSnapshot::generatePosition(int m) const
{
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    _cellv[m].extent(rmin, thetamin, phimin, rmax, thetamax, phimax);

    double r = cbrt(Cubic::pow3(rmin) + Cubic::pow3(rmin, rmax) * random()->uniform());
    double theta = acos(cos(thetamin) + (cos(thetamax) - cos(thetamin)) * random()->uniform());
    double phi = phimin + (phimax - phimin) * random()->uniform();
    return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
}

////////////////////////////////////////////////////////////////////

Position SphericalCellSnapshot::generatePosition() const
{
    // if there are no cells, return the origin
    if (_propv.empty()) return Position();

    // select a cell according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////

const Array& SphericalCellSnapshot::properties(int m) const
{
    return _propv[m];
}

////////////////////////////////////////////////////////////////////

void SphericalCellSnapshot::getEntities(EntityCollection& entities, Position bfr) const
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

void SphericalCellSnapshot::getEntities(EntityCollection& entities, Position bfr, Direction bfk) const
{
    entities.clear();
    for (int m : _search.entitiesFor(bfr, bfk))
    {
        entities.add(m, _cellv[m].intersection(bfr, bfk));
    };
}

////////////////////////////////////////////////////////////////////

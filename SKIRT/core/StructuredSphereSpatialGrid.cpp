/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StructuredSphereSpatialGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Quadrics.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void StructuredSphereSpatialGrid::setupSelfAfter()
{
    SphereSpatialGrid::setupSelfAfter();

    // ---- radial ----

    // set up the radial grid
    _Nr = _meshRadial->numBins();
    _rv = _meshRadial->mesh() * (maxRadius() - minRadius()) + minRadius();

    // limit the epsilon we use for progressing the path to a value smaller than the hole and/or the first bin
    bool hasHole = _rv[0] > 0.;
    _eps = min(meshEpsilonScale() * maxRadius(), 0.1 * (hasHole ? min(_rv[0], _rv[1] - _rv[0]) : _rv[1]));

    // if the grid has no hole, create an artificial hole larger than the epsilon we use for progressing the path
    // so that the segment generator has a chance to reset the phi bin index when crossing the origin
    if (!hasHole) _rv[0] = 2. * _eps;

    // ---- polar ----

    // set up the polar grid; pre-calculate the cosines for each angular boundary
    _Ntheta = initPolarGrid(_thetav, _cv, _meshPolar);

    // ---- azimuthal ----

    // set up the azimuthal grid
    _Nphi = _meshAzimuthal->numBins();
    _phiv = _meshAzimuthal->mesh() * (2. * M_PI) - M_PI;

    // verify that the azimuth bins are smaller than 120 degrees, with some leeway,
    // so that the segment generator never inadvertently intersects the path with the reflected phi bin border
    constexpr double maxPhi = 0.7 * M_PI;
    for (int k = 0; k != _Nphi; ++k)
        if (_phiv[k + 1] - _phiv[k] > maxPhi) throw FATALERROR("Azimuth bin is wider than 120 deg");

    // pre-calculate sines and cosines for azimuthal bin borders; make sure that the outer boundary values are exact
    _sinv = sin(_phiv);
    _cosv = cos(_phiv);
    _sinv[0] = 0.;
    _cosv[0] = -1.;
    _sinv[_Nphi] = 0.;
    _cosv[_Nphi] = -1.;

    // ---- cells ----

    // cache the final number of cells
    _Ncells = _Nr * _Ntheta * _Nphi;

    // precompute the nominal volume of each structured cell
    _cellVolume.resize(_Ncells);
    for (int m = 0; m != _Ncells; ++m)
    {
        double rmin, thetamin, phimin, rmax, thetamax, phimax;
        if (getCoords(m, rmin, thetamin, phimin, rmax, thetamax, phimax))
            _cellVolume[m] =
                (1. / 3.) * Quadrics::cube(rmin, rmax) * (cos(thetamin) - cos(thetamax)) * (phimax - phimin);
    }
}

//////////////////////////////////////////////////////////////////////

double StructuredSphereSpatialGrid::meshEpsilonScale() const
{
    return 1e-12;
}

//////////////////////////////////////////////////////////////////////

int StructuredSphereSpatialGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int StructuredSphereSpatialGrid::numCells() const
{
    return _Ncells;
}

//////////////////////////////////////////////////////////////////////

double StructuredSphereSpatialGrid::volume(int m) const
{
    if (m >= 0 && m < _Ncells) return _cellVolume[m];
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double StructuredSphereSpatialGrid::diagonal(int m) const
{
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(m, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        Position p0(rmin, thetamin, phimin, Position::CoordinateSystem::SPHERICAL);
        Position p1(rmax, thetamax, phimax, Position::CoordinateSystem::SPHERICAL);
        return (p1 - p0).norm();
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

int StructuredSphereSpatialGrid::cellIndex(Position bfr) const
{
    double r, theta, phi;
    bfr.spherical(r, theta, phi);

    int i = NR::locateFail(_rv, r);
    if (i < 0) return -1;
    int j = NR::locateClip(_thetav, theta);
    int k = NR::locateClip(_phiv, phi);

    return index(i, j, k);
}

//////////////////////////////////////////////////////////////////////

Position StructuredSphereSpatialGrid::centralPositionInCell(int m) const
{
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(m, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        double r = 0.5 * (rmin + rmax);
        double theta = 0.5 * (thetamin + thetamax);
        double phi = 0.5 * (phimin + phimax);
        return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

Position StructuredSphereSpatialGrid::randomPositionInCell(int m) const
{
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(m, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        double r = cbrt(Quadrics::cube(rmin) + Quadrics::cube(rmin, rmax) * random()->uniform());
        double theta = acos(cos(thetamin) + (cos(thetamax) - cos(thetamin)) * random()->uniform());
        double phi = phimin + (phimax - phimin) * random()->uniform();
        return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

void StructuredSphereSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // spheres
    for (double r : _rv) outfile->writeCircle(r);

    // meridional planes
    for (double phi : _phiv)
        outfile->writeLine(_rv[0] * cos(phi), _rv[0] * sin(phi), _rv[_Nr] * cos(phi), _rv[_Nr] * sin(phi));
}

//////////////////////////////////////////////////////////////////////

void StructuredSphereSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    writeMeridionalStructure(outfile);
}

//////////////////////////////////////////////////////////////////////

void StructuredSphereSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    writeMeridionalStructure(outfile);
}

//////////////////////////////////////////////////////////////////////

void StructuredSphereSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // for all spheres
    for (double r : _rv)
    {
        // draw the intersections of the spheres with the cones
        for (double theta : _thetav) outfile->writeCircle(r * sin(theta), r * cos(theta));

        // draw the intersections of the spheres with the meridional planes
        for (double phi : _phiv) outfile->writeMeridionalHalfCircle(r, phi);
    }

    // draw the intersections of the cones with the meridional planes
    for (double theta : _thetav)
    {
        for (double phi : _phiv)
        {
            outfile->writeLine(_rv[0] * sin(theta) * cos(phi), _rv[0] * sin(theta) * sin(phi), _rv[0] * cos(theta),
                               _rv[_Nr] * sin(theta) * cos(phi), _rv[_Nr] * sin(theta) * sin(phi),
                               _rv[_Nr] * cos(theta));
        }
    }
}

//////////////////////////////////////////////////////////////////////

void StructuredSphereSpatialGrid::writeMeridionalStructure(SpatialGridPlotFile* outfile) const
{
    // spheres
    for (double r : _rv) outfile->writeCircle(r);

    // cones
    for (double theta : _thetav)
    {
        outfile->writeLine(_rv[0] * sin(theta), _rv[0] * cos(theta), _rv[_Nr] * sin(theta), _rv[_Nr] * cos(theta));
        outfile->writeLine(-_rv[0] * sin(theta), -_rv[0] * cos(theta), -_rv[_Nr] * sin(theta), -_rv[_Nr] * cos(theta));
    }
}

//////////////////////////////////////////////////////////////////////

int StructuredSphereSpatialGrid::index(int i, int j, int k) const
{
    return k + (j + i * _Ntheta) * _Nphi;
}

//////////////////////////////////////////////////////////////////////

bool StructuredSphereSpatialGrid::getCoords(int m, double& rmin, double& thetamin, double& phimin, double& rmax,
                                            double& thetamax, double& phimax) const
{
    if (m < 0 || m >= _Ncells) return false;

    int i = m / (_Ntheta * _Nphi);
    int j = (m / _Nphi) % _Ntheta;
    int k = m % _Nphi;

    rmin = _rv[i];
    thetamin = _thetav[j];
    phimin = _phiv[k];
    rmax = _rv[i + 1];
    thetamax = _thetav[j + 1];
    phimax = _phiv[k + 1];
    return true;
}

//////////////////////////////////////////////////////////////////////

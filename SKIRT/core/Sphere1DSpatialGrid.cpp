/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Sphere1DSpatialGrid.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void Sphere1DSpatialGrid::setupSelfAfter()
{
    // Set up the grid properties
    _Nr = _meshRadial->numBins();
    _rv = minRadius() + _meshRadial->mesh() * (maxRadius() - minRadius());

    // base class setupSelfAfter() depends on initialization performed above
    SphereSpatialGrid::setupSelfAfter();
}

//////////////////////////////////////////////////////////////////////

int Sphere1DSpatialGrid::dimension() const
{
    return 1;
}

//////////////////////////////////////////////////////////////////////

int Sphere1DSpatialGrid::numCells() const
{
    return _Nr;
}

//////////////////////////////////////////////////////////////////////

double Sphere1DSpatialGrid::volume(int m) const
{
    int i = m;
    if (i < 0 || i >= _Nr)
        return 0.0;
    else
    {
        double rL = _rv[i];
        double rR = _rv[i + 1];
        return 4.0 * M_PI / 3.0 * (rR - rL) * (rR * rR + rR * rL + rL * rL);
    }
}

//////////////////////////////////////////////////////////////////////

int Sphere1DSpatialGrid::cellIndex(Position bfr) const
{
    return NR::locateFail(_rv, bfr.radius());
}

//////////////////////////////////////////////////////////////////////

Position Sphere1DSpatialGrid::centralPositionInCell(int m) const
{
    int i = m;
    double r = (_rv[i] + _rv[i + 1]) / 2.0;
    return Position(r, 0, 0);
}

//////////////////////////////////////////////////////////////////////

Position Sphere1DSpatialGrid::randomPositionInCell(int m) const
{
    int i = m;
    Direction bfk = random()->direction();
    double r = _rv[i] + (_rv[i + 1] - _rv[i]) * random()->uniform();
    return Position(r, bfk);
}

//////////////////////////////////////////////////////////////////////

void Sphere1DSpatialGrid::path(SpatialGridPath* path) const
{
    // Determination of the initial position and direction of the path,
    // and calculation of some initial values

    path->clear();
    double x, y, z;
    path->position().cartesian(x, y, z);
    double kx, ky, kz;
    path->direction().cartesian(kx, ky, kz);
    double rmax = maxRadius();
    double rmin = minRadius();

    // Move the photon packet to the first grid cell that it will pass.
    // If it does not pass any grid cell, return an empty path.

    double r = path->position().radius();
    double q = x * kx + y * ky + z * kz;
    double p = sqrt((r - q) * (r + q));
    if (r > rmax)
    {
        if (q > 0.0 || p > rmax)
            return path->clear();
        else
        {
            r = rmax - 1e-8 * (_rv[_Nr] - _rv[_Nr - 1]);
            // path intersects rmax boundary going inward, qmax is negative
            double qmax = -sqrt((rmax - p) * (rmax + p));
            double ds = (qmax - q);
            path->addSegment(-1, ds);
            q = qmax;
        }
    }
    else if (r < rmin)
    {
        r = rmin + 1e-8 * (_rv[1] - _rv[0]);
        // path intersects rmin boundary going outward, qmin is positive
        double qmin = sqrt((rmin - p) * (rmin + p));
        double ds = (qmin - q);
        path->addSegment(-1, ds);
        q = qmin;
    }

    // Determination of the initial grid cell

    int i = NR::locateClip(_rv, r);

    // And here we go...

    double rN, qN;

    // Inward movement (not always...)

    if (q < 0.0)
    {
        int imin = NR::locate(_rv, p);  // returns -1 when p < rmin
        while (i > imin)
        {
            rN = _rv[i];  // i >= 0 here
            qN = -sqrt((rN - p) * (rN + p));
            int m = i;
            double ds = qN - q;
            path->addSegment(m, ds);
            i--;  // i can become -1 here
            q = qN;
        }
    }

    // Outward movement (starting with i potentially -1)

    while (true)
    {
        rN = _rv[i + 1];
        qN = sqrt((rN - p) * (rN + p));
        int m = i;
        double ds = qN - q;
        path->addSegment(m, ds);  // m can be -1 here, as intended
        i++;
        if (i >= _Nr)
            return;
        else
            q = qN;
    }
}

//////////////////////////////////////////////////////////////////////

void Sphere1DSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _Nr; i++) outfile->writeCircle(_rv[i]);
}

//////////////////////////////////////////////////////////////////////

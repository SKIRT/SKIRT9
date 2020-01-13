/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Cylinder2DSpatialGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void Cylinder2DSpatialGrid::setupSelfAfter()
{
    // initialize our local mesh arrays
    _NR = _meshRadial->numBins();
    _Nz = _meshZ->numBins();
    double Rmax = maxRadius();
    double zmin = minZ();
    double zmax = maxZ();
    _Rv = _meshRadial->mesh() * Rmax;
    _zv = _meshZ->mesh() * (zmax - zmin) + zmin;

    // base class setupSelfAfter() depends on initialization performed above
    CylinderSpatialGrid::setupSelfAfter();
}

//////////////////////////////////////////////////////////////////////

int Cylinder2DSpatialGrid::dimension() const
{
    return 2;
}

//////////////////////////////////////////////////////////////////////

int Cylinder2DSpatialGrid::numCells() const
{
    return _NR * _Nz;
}

//////////////////////////////////////////////////////////////////////

double Cylinder2DSpatialGrid::volume(int m) const
{
    int i, k;
    invertIndex(m, i, k);
    if (i < 0 || i >= _NR || k < 0 || k >= _Nz)
        return 0.0;
    else
        return M_PI * (_zv[k + 1] - _zv[k]) * (_Rv[i + 1] - _Rv[i]) * (_Rv[i + 1] + _Rv[i]);
}

//////////////////////////////////////////////////////////////////////

int Cylinder2DSpatialGrid::cellIndex(Position bfr) const
{
    int i = NR::locateFail(_Rv, bfr.cylRadius());
    int k = NR::locateFail(_zv, bfr.height());
    if (i < 0 || k < 0) return -1;
    return index(i, k);
}

//////////////////////////////////////////////////////////////////////

Position Cylinder2DSpatialGrid::centralPositionInCell(int m) const
{
    int i, k;
    invertIndex(m, i, k);
    double R = (_Rv[i] + _Rv[i + 1]) / 2.0;
    double phi = 0.0;
    double z = (_zv[k] + _zv[k + 1]) / 2.0;
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////

Position Cylinder2DSpatialGrid::randomPositionInCell(int m) const
{
    int i, k;
    invertIndex(m, i, k);
    double R = _Rv[i] + (_Rv[i + 1] - _Rv[i]) * random()->uniform();
    double phi = 2.0 * M_PI * random()->uniform();
    double z = _zv[k] + (_zv[k + 1] - _zv[k]) * random()->uniform();
    return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
}

//////////////////////////////////////////////////////////////////////

void Cylinder2DSpatialGrid::path(SpatialGridPath* path) const
{
    // Determination of the initial position and direction of the path,
    // and calculation of some initial values
    path->clear();
    double kx, ky, kz;
    path->direction().cartesian(kx, ky, kz);
    double kq = sqrt(kx * kx + ky * ky);
    if (kz == 0.0) kz = 1e-20;  // avoid moving exactly parallel to the equatorial plane
    if (kq == 0.0) kq = 1e-20;  // avoid moving exactly parallel to the z-axis
    double x, y, z;
    path->position().cartesian(x, y, z);
    double R = path->position().cylRadius();
    double q = (x * kx + y * ky) / kq;
    double p2 = (R - q) * (R + q);
    double p = sqrt(max(0.0, p2));  // make sure that p>=0 here; necessary sometimes due to rounding errors
    double Rmax = maxRadius();
    double zmin = minZ();
    double zmax = maxZ();

    // Move the photon packet to the first grid cell that it will pass.
    // If it does not pass any grid cell, return an empty path.
    // Otherwise calculate the distance covered and add a segment to the path.

    if (R >= Rmax)
    {
        if (q > 0.0 || p > Rmax)
            return path->clear();
        else
        {
            R = Rmax - 1e-8 * (_Rv[_NR] - _Rv[_NR - 1]);
            double qmax = sqrt((Rmax - p) * (Rmax + p));
            double ds = (qmax - q) / kq;
            path->addSegment(-1, ds);
            q = qmax;
            z += kz * ds;
        }
    }
    if (z < zmin)
    {
        if (kz <= 0.0)
            return path->clear();
        else
        {
            double ds = (zmin - z) / kz;
            path->addSegment(-1, ds);
            q += kq * ds;
            R = sqrt(p * p + q * q);
            z = zmin + 1e-8 * (_zv[1] - _zv[0]);
        }
    }
    else if (z > zmax)
    {
        if (kz >= 0.0)
            return path->clear();
        else
        {
            double ds = (zmax - z) / kz;
            path->addSegment(-1, ds);
            q += kq * ds;
            R = sqrt(p * p + q * q);
            z = zmax - 1e-8 * (_zv[_Nz] - _zv[_Nz - 1]);
        }
    }
    if (std::isinf(R) || std::isnan(R) || std::isinf(z) || std::isnan(z) || R >= Rmax || z <= zmin || z >= zmax)
        return path->clear();

    // Determination of the initial grid cell

    int i = NR::locateClip(_Rv, R);
    int k = NR::locateClip(_zv, z);

    // And here we go...

    double ds, dsq, dsz;
    double RN, qN, zN;

    // SCENARIO 1: UPWARD MOVEMENT

    if (kz >= 0.0)
    {
        if (q < 0.0)
        {
            int imin = NR::locateClip(_Rv, p);
            RN = _Rv[i];
            qN = -sqrt((RN - p) * (RN + p));
            zN = _zv[k + 1];
            while (i > imin)
            {
                int m = index(i, k);
                dsq = (qN - q) / kq;
                dsz = (zN - z) / kz;
                if (dsq < dsz)
                {
                    ds = dsq;
                    path->addSegment(m, ds);
                    i--;
                    q = qN;
                    z += kz * ds;
                    RN = _Rv[i];
                    qN = -sqrt((RN - p) * (RN + p));
                }
                else
                {
                    ds = dsz;
                    path->addSegment(m, ds);
                    k++;
                    if (k >= _Nz)
                        return;
                    else
                    {
                        q += kq * ds;
                        z = zN;
                        zN = _zv[k + 1];
                    }
                }
            }
        }
        RN = _Rv[i + 1];
        qN = sqrt((RN - p) * (RN + p));
        zN = _zv[k + 1];
        while (true)
        {
            int m = index(i, k);
            dsq = (qN - q) / kq;
            dsz = (zN - z) / kz;
            if (dsq < dsz)
            {
                ds = dsq;
                path->addSegment(m, ds);
                i++;
                if (i >= _NR)
                    return;
                else
                {
                    q = qN;
                    z += kz * ds;
                    RN = _Rv[i + 1];
                    qN = sqrt((RN - p) * (RN + p));
                }
            }
            else
            {
                ds = dsz;
                path->addSegment(m, ds);
                k++;
                if (k >= _Nz)
                    return;
                else
                {
                    q += kq * ds;
                    z = zN;
                    zN = _zv[k + 1];
                }
            }
        }
    }

    // SCENARIO 2: DOWNWARD MOVEMENT

    else
    {
        if (q < 0.0)
        {
            int imin = NR::locateClip(_Rv, p);
            RN = _Rv[i];
            qN = -sqrt((RN - p) * (RN + p));
            zN = _zv[k];
            while (i > imin)
            {
                int m = index(i, k);
                dsq = (qN - q) / kq;
                dsz = (zN - z) / kz;
                if (dsq < dsz)
                {
                    ds = dsq;
                    path->addSegment(m, ds);
                    i--;
                    q = qN;
                    z += kz * ds;
                    RN = _Rv[i];
                    qN = -sqrt((RN - p) * (RN + p));
                }
                else
                {
                    ds = dsz;
                    path->addSegment(m, ds);
                    k--;
                    if (k < 0)
                        return;
                    else
                    {
                        q += kq * ds;
                        z = zN;
                        zN = _zv[k];
                    }
                }
            }
        }
        RN = _Rv[i + 1];
        qN = sqrt((RN - p) * (RN + p));
        zN = _zv[k];
        while (true)
        {
            int m = index(i, k);
            dsq = (qN - q) / kq;
            dsz = (zN - z) / kz;
            if (dsq < dsz)
            {
                ds = dsq;
                path->addSegment(m, ds);
                i++;
                if (i >= _NR)
                    return;
                else
                {
                    q = qN;
                    z += kz * ds;
                    RN = _Rv[i + 1];
                    qN = sqrt((RN - p) * (RN + p));
                }
            }
            else
            {
                ds = dsz;
                path->addSegment(m, ds);
                k--;
                if (k < 0)
                    return;
                else
                {
                    q += kq * ds;
                    z = zN;
                    zN = _zv[k];
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////

void Cylinder2DSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _NR; i++) outfile->writeCircle(_Rv[i]);
}

//////////////////////////////////////////////////////////////////////

void Cylinder2DSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    double Rmax = maxRadius();
    double zmin = minZ();
    double zmax = maxZ();
    for (int i = 0; i <= _NR; i++)
    {
        outfile->writeLine(_Rv[i], zmin, _Rv[i], zmax);
        outfile->writeLine(-_Rv[i], zmin, -_Rv[i], zmax);
    }
    for (int k = 0; k <= _Nz; k++)
    {
        outfile->writeLine(-Rmax, _zv[k], Rmax, _zv[k]);
    }
}

//////////////////////////////////////////////////////////////////////

int Cylinder2DSpatialGrid::index(int i, int k) const
{
    return k + _Nz * i;
}

//////////////////////////////////////////////////////////////////////

void Cylinder2DSpatialGrid::invertIndex(int m, int& i, int& k) const
{
    i = m / _Nz;
    k = m % _Nz;
}

//////////////////////////////////////////////////////////////////////

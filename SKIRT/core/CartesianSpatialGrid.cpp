/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CartesianSpatialGrid.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::setupSelfAfter()
{
    // initialize our local mesh arrays
    _Nx = _meshX->numBins();
    _Ny = _meshY->numBins();
    _Nz = _meshZ->numBins();
    _xv = _meshX->mesh()*(xmax()-xmin()) + xmin();
    _yv = _meshY->mesh()*(ymax()-ymin()) + ymin();
    _zv = _meshZ->mesh()*(zmax()-zmin()) + zmin();

    // base class setupSelfAfter() depends on initialization performed above
    BoxSpatialGrid::setupSelfAfter();
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::numCells() const
{
    return _Nx*_Ny*_Nz;
}

//////////////////////////////////////////////////////////////////////

double CartesianSpatialGrid::volume(int m) const
{
    return box(m).volume();
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::cellIndex(Position bfr) const
{
    int i = NR::locateFail(_xv, bfr.x());
    int j = NR::locateFail(_yv, bfr.y());
    int k = NR::locateFail(_zv, bfr.z());
    if (i<0 || j<0 || k<0)
        return -1;
    else
        return index(i,j,k);
}

//////////////////////////////////////////////////////////////////////

Position CartesianSpatialGrid::centralPositionInCell(int m) const
{
    return Position(box(m).center());
}

//////////////////////////////////////////////////////////////////////

Position CartesianSpatialGrid::randomPositionInCell(int m) const
{
    return random()->position(box(m));
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::path(SpatialGridPath* path) const
{
    // Determination of the initial position and direction of the path,
    // and calculation of some initial values

    path->clear();
    double kx,ky,kz;
    path->direction().cartesian(kx,ky,kz);
    double x,y,z;
    path->position().cartesian(x,y,z);
    double ds, dsx, dsy, dsz;

    // move the photon package to the first grid cell that it will
    // pass. If it does not pass any grid cell, return an empty path.

    if (x<xmin())
    {
        if (kx<=0.0) return path->clear();
        else
        {
            ds = (xmin()-x)/kx;
            path->addSegment(-1,ds);
            x = xmin() + 1e-8*(_xv[1]-_xv[0]);
            y += ky*ds;
            z += kz*ds;
        }
    }
    else if (x>xmax())
    {
        if (kx>=0.0) return path->clear();
        else
        {
            ds = (xmax()-x)/kx;
            path->addSegment(-1,ds);
            x = xmax() - 1e-8*(_xv[_Nx]-_xv[_Nx-1]);
            y += ky*ds;
            z += kz*ds;
        }
    }
    if (y<ymin())
    {
        if (ky<=0.0) return path->clear();
        else
        {
            ds = (ymin()-y)/ky;
            path->addSegment(-1,ds);
            x += kx*ds;
            y = ymin() + 1e-8*(_yv[1]-_yv[0]);
            z += kz*ds;
        }
    }
    else if (y>ymax())
    {
        if (ky>=0.0) return path->clear();
        else
        {
            ds = (ymax()-y)/ky;
            path->addSegment(-1,ds);
            x += kx*ds;
            y = ymax() - 1e-8*(_yv[_Ny]-_yv[_Ny-1]);
            z += kz*ds;
        }
    }
    if (z<zmin())
    {
        if (kz<=0.0) return path->clear();
        else
        {
            ds = (zmin()-z)/kz;
            path->addSegment(-1,ds);
            x += kx*ds;
            y += ky*ds;
            z = zmin() + 1e-8*(_zv[1]-_zv[0]);
        }
    }
    else if (z>zmax())
    {
        if (kz>=0.0) return path->clear();
        else
        {
            ds = (zmax()-z)/kz;
            path->addSegment(-1,ds);
            x += kx*ds;
            y += ky*ds;
            z = zmax() - 1e-8*(_zv[_Nz]-_zv[_Nz-1]);
        }
    }

    if (!contains(x,y,z)) return path->clear();

    // Now determine which grid cell we are in...

    int i = NR::locateClip(_xv, x);
    int j = NR::locateClip(_yv, y);
    int k = NR::locateClip(_zv, z);

    // And there we go...

    while (true)
    {
        int m = index(i,j,k);
        double xE = (kx<0.0) ? _xv[i] : _xv[i+1];
        double yE = (ky<0.0) ? _yv[j] : _yv[j+1];
        double zE = (kz<0.0) ? _zv[k] : _zv[k+1];
        dsx = (fabs(kx)>1e-15) ? (xE-x)/kx : DBL_MAX;
        dsy = (fabs(ky)>1e-15) ? (yE-y)/ky : DBL_MAX;
        dsz = (fabs(kz)>1e-15) ? (zE-z)/kz : DBL_MAX;
        if (dsx<=dsy && dsx<=dsz)
        {
            ds = dsx;
            path->addSegment(m, ds);
            i += (kx<0.0) ? -1 : 1;
            if (i>=_Nx || i<0) return;
            else
            {
                x = xE;
                y += ky*ds;
                z += kz*ds;
            }
        }
        else if (dsy< dsx && dsy<=dsz)
        {
            ds = dsy;
            path->addSegment(m, ds);
            j += (ky<0.0) ? -1 : 1;
            if (j>=_Ny || j<0) return;
            else
            {
                x += kx*ds;
                y = yE;
                z += kz*ds;
            }
        }
        else if (dsz< dsx && dsz< dsy)
        {
            ds = dsz;
            path->addSegment(m, ds);
            k += (kz<0.0) ? -1 : 1;
            if (k>=_Nz || k<0) return;
            else
            {
                x += kx*ds;
                y += ky*ds;
                z  = zE;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    for (int i=0; i<=_Nx; i++) outfile->writeLine(_xv[i], ymin(), _xv[i], ymax());
    for (int j=0; j<=_Ny; j++) outfile->writeLine(xmin(), _yv[j], xmax(), _yv[j]);
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    for (int i=0; i<=_Nx; i++) outfile->writeLine(_xv[i], zmin(), _xv[i], zmax());
    for (int k=0; k<=_Nz; k++) outfile->writeLine(xmin(), _zv[k], xmax(), _zv[k]);
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    for (int j=0; j<=_Ny; j++) outfile->writeLine(_yv[j], zmin(), _yv[j], zmax());
    for (int k=0; k<=_Nz; k++) outfile->writeLine(ymin(), _zv[k], ymax(), _zv[k]);
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    for (int i=0; i<=_Nx; i++)
        for (int j=0; j<=_Ny; j++)
            outfile->writeLine(_xv[i], _yv[j], zmin(), _xv[i], _yv[j], zmax());
    for (int i=0; i<=_Nx; i++)
        for (int k=0; k<=_Nz; k++)
            outfile->writeLine(_xv[i], ymin(), _zv[k], _xv[i], ymax(), _zv[k]);
    for (int j=0; j<=_Ny; j++)
        for (int k=0; k<=_Nz; k++)
            outfile->writeLine(xmin(), _yv[j], _zv[k], xmax(), _yv[j], _zv[k]);
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::index(int i, int j, int k) const
{
    return k+_Nz*j+_Nz*_Ny*i;
}

//////////////////////////////////////////////////////////////////////

Box CartesianSpatialGrid::box(int m) const
{
    int i = m / (_Nz*_Ny);
    int j = (m/_Nz) % _Ny;
    int k = m % _Nz;

    if (i<0 || j<0 || k<0 || i>=_Nx || j>=_Ny || k>=_Nz)
        return Box();
    else
        return Box(_xv[i],_yv[j],_zv[k],_xv[i+1],_yv[j+1],_zv[k+1]);
}

//////////////////////////////////////////////////////////////////////

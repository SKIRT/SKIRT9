/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CartesianSpatialGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"
#include "SpatialGridSegmentGenerator.hpp"

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::setupSelfAfter()
{
    // initialize our local mesh arrays
    _Nx = _meshX->numBins();
    _Ny = _meshY->numBins();
    _Nz = _meshZ->numBins();
    _xv = _meshX->mesh() * (xmax() - xmin()) + xmin();
    _yv = _meshY->mesh() * (ymax() - ymin()) + ymin();
    _zv = _meshZ->mesh() * (zmax() - zmin()) + zmin();

    // base class setupSelfAfter() depends on initialization performed above
    BoxSpatialGrid::setupSelfAfter();
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::numCells() const
{
    return _Nx * _Ny * _Nz;
}

//////////////////////////////////////////////////////////////////////

double CartesianSpatialGrid::volume(int m) const
{
    return box(m).volume();
}

//////////////////////////////////////////////////////////////////////

double CartesianSpatialGrid::diagonal(int m) const
{
    return box(m).diagonal();
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::cellIndex(Position bfr) const
{
    int i = NR::locateFail(_xv, bfr.x());
    int j = NR::locateFail(_yv, bfr.y());
    int k = NR::locateFail(_zv, bfr.z());
    if (i < 0 || j < 0 || k < 0)
        return -1;
    else
        return index(i, j, k);
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

class CartesianSpatialGrid::SegmentGenerator : public SpatialGridSegmentGenerator
{
    const CartesianSpatialGrid* _grid{nullptr};
    int _i{-1}, _j{-1}, _k{-1};
    enum class CurrentPosition { Unknown, Inside, Outside };
    CurrentPosition _curpos{CurrentPosition::Unknown};

public:
    SegmentGenerator(const SpatialGridPath* path, const CartesianSpatialGrid* grid)
        : SpatialGridSegmentGenerator(path), _grid(grid)
    {}

    bool next() override
    {
        switch (_curpos)
        {
            case CurrentPosition::Unknown:
            {
                // move the photon packet inside the grid; if this is impossible, return an empty path
                moveInside(_grid->extent(), 1e-12 * _grid->extent().widths().norm());
                if (!_grid->contains(r())) return false;

                // determine which grid cell we are in
                _i = NR::locateClip(_grid->_xv, rx());
                _j = NR::locateClip(_grid->_yv, ry());
                _k = NR::locateClip(_grid->_zv, rz());
                _curpos = CurrentPosition::Inside;

                // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
                // otherwise fall through to determine the first actual segment
                if (ds() > 0.) return true;
            }

            case CurrentPosition::Inside:
            {
                // determine the segment from the current position to the first cell wall
                // and adjust the position and cell indices accordingly
                int m = _grid->index(_i, _j, _k);
                double xE = (kx() < 0.0) ? _grid->_xv[_i] : _grid->_xv[_i + 1];
                double yE = (ky() < 0.0) ? _grid->_yv[_j] : _grid->_yv[_j + 1];
                double zE = (kz() < 0.0) ? _grid->_zv[_k] : _grid->_zv[_k + 1];
                double dsx = (fabs(kx()) > 1e-15) ? (xE - rx()) / kx() : DBL_MAX;
                double dsy = (fabs(ky()) > 1e-15) ? (yE - ry()) / ky() : DBL_MAX;
                double dsz = (fabs(kz()) > 1e-15) ? (zE - rz()) / kz() : DBL_MAX;

                if (dsx <= dsy && dsx <= dsz)
                {
                    setSegment(m, dsx);
                    setrx(xE);
                    propagatey();
                    propagatez();
                    _i += (kx() < 0.0) ? -1 : 1;
                    if (_i >= _grid->_Nx || _i < 0) _curpos = CurrentPosition::Outside;
                }
                else if (dsy < dsx && dsy <= dsz)
                {
                    setSegment(m, dsy);
                    setry(yE);
                    propagatex();
                    propagatez();
                    _j += (ky() < 0.0) ? -1 : 1;
                    if (_j >= _grid->_Ny || _j < 0) _curpos = CurrentPosition::Outside;
                }
                else if (dsz < dsx && dsz < dsy)
                {
                    setSegment(m, dsz);
                    setrz(zE);
                    propagatex();
                    propagatey();
                    _k += (kz() < 0.0) ? -1 : 1;
                    if (_k >= _grid->_Nz || _k < 0) _curpos = CurrentPosition::Outside;
                }
                else
                {
                    throw FATALERROR("path segment error");
                }
                return true;
            }

            case CurrentPosition::Outside:
            {
            }
        }
        return false;
    }
};

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::path(SpatialGridPath* path) const
{
    CartesianSpatialGrid::SegmentGenerator generator(path, this);
    path->clear();
    while (generator.next())
    {
        path->addSegment(generator.m(), generator.ds());
    }

/*
    // If the photon packet starts outside the grid, move it inside;
    // if this is impossible, return an empty path
    double eps = 1e-12 * extent().widths().norm();
    Position bfr = path->moveInside(extent(), eps);
    if (!contains(bfr)) return path->clear();

    // Get the direction and the current position of the path
    double kx, ky, kz;
    path->direction().cartesian(kx, ky, kz);
    double x, y, z;
    bfr.cartesian(x, y, z);

    // Determine which grid cell we are in
    int i = NR::locateClip(_xv, x);
    int j = NR::locateClip(_yv, y);
    int k = NR::locateClip(_zv, z);

    // There we go...
    double ds, dsx, dsy, dsz;
    while (true)
    {
        int m = index(i, j, k);
        double xE = (kx < 0.0) ? _xv[i] : _xv[i + 1];
        double yE = (ky < 0.0) ? _yv[j] : _yv[j + 1];
        double zE = (kz < 0.0) ? _zv[k] : _zv[k + 1];
        dsx = (fabs(kx) > 1e-15) ? (xE - x) / kx : DBL_MAX;
        dsy = (fabs(ky) > 1e-15) ? (yE - y) / ky : DBL_MAX;
        dsz = (fabs(kz) > 1e-15) ? (zE - z) / kz : DBL_MAX;
        if (dsx <= dsy && dsx <= dsz)
        {
            ds = dsx;
            path->addSegment(m, ds);
            i += (kx < 0.0) ? -1 : 1;
            if (i >= _Nx || i < 0)
                return;
            else
            {
                x = xE;
                y += ky * ds;
                z += kz * ds;
            }
        }
        else if (dsy < dsx && dsy <= dsz)
        {
            ds = dsy;
            path->addSegment(m, ds);
            j += (ky < 0.0) ? -1 : 1;
            if (j >= _Ny || j < 0)
                return;
            else
            {
                x += kx * ds;
                y = yE;
                z += kz * ds;
            }
        }
        else if (dsz < dsx && dsz < dsy)
        {
            ds = dsz;
            path->addSegment(m, ds);
            k += (kz < 0.0) ? -1 : 1;
            if (k >= _Nz || k < 0)
                return;
            else
            {
                x += kx * ds;
                y += ky * ds;
                z = zE;
            }
        }
    }
*/
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _Nx; i++) outfile->writeLine(_xv[i], ymin(), _xv[i], ymax());
    for (int j = 0; j <= _Ny; j++) outfile->writeLine(xmin(), _yv[j], xmax(), _yv[j]);
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _Nx; i++) outfile->writeLine(_xv[i], zmin(), _xv[i], zmax());
    for (int k = 0; k <= _Nz; k++) outfile->writeLine(xmin(), _zv[k], xmax(), _zv[k]);
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    for (int j = 0; j <= _Ny; j++) outfile->writeLine(_yv[j], zmin(), _yv[j], zmax());
    for (int k = 0; k <= _Nz; k++) outfile->writeLine(ymin(), _zv[k], ymax(), _zv[k]);
}

//////////////////////////////////////////////////////////////////////

void CartesianSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _Nx; i++)
        for (int j = 0; j <= _Ny; j++) outfile->writeLine(_xv[i], _yv[j], zmin(), _xv[i], _yv[j], zmax());
    for (int i = 0; i <= _Nx; i++)
        for (int k = 0; k <= _Nz; k++) outfile->writeLine(_xv[i], ymin(), _zv[k], _xv[i], ymax(), _zv[k]);
    for (int j = 0; j <= _Ny; j++)
        for (int k = 0; k <= _Nz; k++) outfile->writeLine(xmin(), _yv[j], _zv[k], xmax(), _yv[j], _zv[k]);
}

//////////////////////////////////////////////////////////////////////

int CartesianSpatialGrid::index(int i, int j, int k) const
{
    return k + _Nz * j + _Nz * _Ny * i;
}

//////////////////////////////////////////////////////////////////////

Box CartesianSpatialGrid::box(int m) const
{
    int i = m / (_Nz * _Ny);
    int j = (m / _Nz) % _Ny;
    int k = m % _Nz;

    if (i < 0 || j < 0 || k < 0 || i >= _Nx || j >= _Ny || k >= _Nz)
        return Box();
    else
        return Box(_xv[i], _yv[j], _zv[k], _xv[i + 1], _yv[j + 1], _zv[k + 1]);
}

//////////////////////////////////////////////////////////////////////

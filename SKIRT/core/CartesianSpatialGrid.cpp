/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CartesianSpatialGrid.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"

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

class CartesianSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const CartesianSpatialGrid* _grid{nullptr};
    int _i{-1}, _j{-1}, _k{-1};

public:
    MySegmentGenerator(const CartesianSpatialGrid* grid) : _grid(grid) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // try moving the photon packet inside the grid; if this is impossible, return an empty path
                if (!moveInside(_grid->extent(), 1e-12 * _grid->extent().widths().norm())) return false;

                // determine which grid cell we are in
                _i = NR::locateClip(_grid->_xv, rx());
                _j = NR::locateClip(_grid->_yv, ry());
                _k = NR::locateClip(_grid->_zv, rz());

                // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
                // otherwise fall through to determine the first actual segment
                if (ds() > 0.) return true;
            }

            // intentionally falls through
            case State::Inside:
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
                    propagatery(dsx);
                    propagaterz(dsx);
                    _i += (kx() < 0.0) ? -1 : 1;
                    if (_i >= _grid->_Nx || _i < 0) setState(State::Outside);
                }
                else if (dsy < dsx && dsy <= dsz)
                {
                    setSegment(m, dsy);
                    setry(yE);
                    propagaterx(dsy);
                    propagaterz(dsy);
                    _j += (ky() < 0.0) ? -1 : 1;
                    if (_j >= _grid->_Ny || _j < 0) setState(State::Outside);
                }
                else  // if (dsz < dsx && dsz < dsy)
                {
                    setSegment(m, dsz);
                    setrz(zE);
                    propagaterx(dsz);
                    propagatery(dsz);
                    _k += (kz() < 0.0) ? -1 : 1;
                    if (_k >= _grid->_Nz || _k < 0) setState(State::Outside);
                }
                return true;
            }

            case State::Outside:
            {
            }
        }
        return false;
    }
};

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> CartesianSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
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

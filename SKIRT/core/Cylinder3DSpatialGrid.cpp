/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Cylinder3DSpatialGrid.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void Cylinder3DSpatialGrid::setupSelfAfter()
{
    CylinderSpatialGrid::setupSelfAfter();

    // initialize our local mesh arrays
    _NR = _meshRadial->numBins();
    _Nphi = _meshAzimuthal->numBins();
    _Nz = _meshZ->numBins();
    _Ncells = _NR * _Nphi * _Nz;
    _Rv = _meshRadial->mesh() * (maxRadius() - minRadius()) + minRadius();
    _phiv = _meshAzimuthal->mesh() * (2. * M_PI) - M_PI;
    _zv = _meshZ->mesh() * (maxZ() - minZ()) + minZ();
}

//////////////////////////////////////////////////////////////////////

int Cylinder3DSpatialGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int Cylinder3DSpatialGrid::numCells() const
{
    return _Ncells;
}

//////////////////////////////////////////////////////////////////////

double Cylinder3DSpatialGrid::volume(int m) const
{
    double Rmin, phimin, zmin, Rmax, phimax, zmax;
    if (getCoords(m, Rmin, phimin, zmin, Rmax, phimax, zmax))
    {
        return 0.5 * (Rmax - Rmin) * (Rmax + Rmin) * (phimax - phimin) * (zmax - zmin);
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double Cylinder3DSpatialGrid::diagonal(int m) const
{
    double Rmin, phimin, zmin, Rmax, phimax, zmax;
    if (getCoords(m, Rmin, phimin, zmin, Rmax, phimax, zmax))
    {
        Position p0(Rmin, phimin, zmin, Position::CoordinateSystem::CYLINDRICAL);
        Position p1(Rmax, phimax, zmax, Position::CoordinateSystem::CYLINDRICAL);
        return (p1 - p0).norm();
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

int Cylinder3DSpatialGrid::cellIndex(Position bfr) const
{
    double R, phi, z;
    bfr.cylindrical(R, phi, z);

    int i = NR::locateFail(_Rv, R);
    if (i < 0) return -1;
    int j = NR::locateClip(_phiv, phi);
    int k = NR::locateFail(_zv, z);
    if (k < 0) return -1;

    return index(i, j, k);
}

//////////////////////////////////////////////////////////////////////

Position Cylinder3DSpatialGrid::centralPositionInCell(int m) const
{
    double Rmin, phimin, zmin, Rmax, phimax, zmax;
    if (getCoords(m, Rmin, phimin, zmin, Rmax, phimax, zmax))
    {
        double R = 0.5 * (Rmin + Rmax);
        double phi = 0.5 * (phimin + phimax);
        double z = 0.5 * (zmin + zmax);
        return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

Position Cylinder3DSpatialGrid::randomPositionInCell(int m) const
{
    double Rmin, phimin, zmin, Rmax, phimax, zmax;
    if (getCoords(m, Rmin, phimin, zmin, Rmax, phimax, zmax))
    {
        double R = sqrt(Rmin * Rmin + (Rmax - Rmin) * (Rmax + Rmin) * random()->uniform());
        double phi = phimin + (phimax - phimin) * random()->uniform();
        double z = zmin + (zmax - zmin) * random()->uniform();
        return Position(R, phi, z, Position::CoordinateSystem::CYLINDRICAL);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

class Cylinder3DSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const Cylinder3DSpatialGrid* _grid{nullptr};
    enum class Phase { UpInwards, UpOutwards, DownInwards, DownOutwards };
    Phase _phase{Phase::UpInwards};
    double _R{0.}, _z{0.}, _kq{0.}, _kz{0.}, _p{0.}, _q{0.};
    int _i{-1}, _imin{-1}, _k{-1};

public:
    MySegmentGenerator(const Cylinder3DSpatialGrid* grid) : _grid(grid) {}

    // This function initializes the position and direction in cylindrical coordinates (R,z,kq,kz,p,q)
    // and attempts to move the position inside the cylinder, setting an empty path segment and the
    // segment generator state accordingly. The function returns true if the position is now inside.
    bool moveInside()
    {
        // intialize position and direction
        _R = sqrt(rx() * rx() + ry() * ry());
        _z = rz();
        _kq = sqrt(kx() * kx() + ky() * ky());
        _kz = kz();
        if (_kq == 0.0) _kq = 1e-20;  // avoid moving exactly parallel to the z-axis
        if (_kz == 0.0) _kz = 1e-20;  // avoid moving exactly parallel to the equatorial plane
        _q = (rx() * kx() + ry() * ky()) / _kq;
        _p = sqrt(max(0., (_R - _q) * (_R + _q)));  // make sure that p>=0 in case of rounding errors

        // get boundaries
        double Rmax = _grid->maxRadius();
        double zmin = _grid->minZ();
        double zmax = _grid->maxZ();

        // initialize to empty segment with zero length
        setEmptySegment();
        setState(State::Outside);

        // keep track of cumulative length of all subsegments
        double cumds = 0.;

        // --> R direction
        if (_R >= Rmax)
        {
            if (_q > 0.0 || _p > Rmax)
                return false;
            else
            {
                double qmin = -sqrt((Rmax - _p) * (Rmax + _p));
                double ds = (qmin - _q) / _kq;
                _q = qmin;
                _R = Rmax - 1e-8 * (_grid->_Rv[_grid->_NR] - _grid->_Rv[_grid->_NR - 1]);
                _z += _kz * ds;
                cumds += ds;
            }
        }

        // --> z direction
        if (_z < zmin)
        {
            if (_kz <= 0.0)
                return false;
            else
            {
                double ds = (zmin - _z) / _kz;
                _q += _kq * ds;
                _R = sqrt(_p * _p + _q * _q);
                _z = zmin + 1e-8 * (_grid->_zv[1] - _grid->_zv[0]);
                cumds += ds;
            }
        }
        else if (_z > zmax)
        {
            if (_kz >= 0.0)
                return false;
            else
            {
                double ds = (zmax - _z) / _kz;
                _q += _kq * ds;
                _R = sqrt(_p * _p + _q * _q);
                _z = zmax - 1e-8 * (_grid->_zv[_grid->_Nz] - _grid->_zv[_grid->_Nz - 1]);
                cumds += ds;
            }
        }

        // in rare border cases, the position can still be outside of the cylinder or have an invalid value
        if (std::isinf(_R) || std::isnan(_R) || std::isinf(_z) || std::isnan(_z)) return false;
        if (_R >= Rmax || _z <= zmin || _z >= zmax) return false;

        // return the empty segment with the cumulative length
        // (which might be zero if the position was inside the box to begin with)
        setEmptySegment(cumds);
        setState(State::Inside);
        return true;
    }

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // try moving the photon packet inside the grid; if this is impossible, return an empty path
                if (!moveInside()) return false;

                // determine the grid cell we are in
                _i = NR::locateClip(_grid->_Rv, _R);
                _k = NR::locateClip(_grid->_zv, _z);

                // determine the initial direction of movement for each coordinate
                if (_kz >= 0.)
                {
                    _phase = Phase::UpOutwards;
                    if (_q < 0.)
                    {
                        _imin = NR::locateClip(_grid->_Rv, _p);
                        if (_i > _imin) _phase = Phase::UpInwards;
                    }
                }
                else
                {
                    _phase = Phase::DownOutwards;
                    if (_q < 0.)
                    {
                        _imin = NR::locateClip(_grid->_Rv, _p);
                        if (_i > _imin) _phase = Phase::DownInwards;
                    }
                }

                // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
                // otherwise fall through to determine the first actual segment
                if (ds() > 0.) return true;
            }

            // intentionally falls through
            case State::Inside:
            {
                switch (_phase)
                {
                    case Phase::UpInwards:
                    {
                        double RN = _grid->_Rv[_i];
                        double qN = -sqrt((RN - _p) * (RN + _p));
                        double zN = _grid->_zv[_k + 1];

                        int m = _grid->index(_i, 0, _k);
                        double dsq = (qN - _q) / _kq;
                        double dsz = (zN - _z) / _kz;
                        if (dsq < dsz)
                        {
                            double ds = dsq;
                            setSegment(m, ds);
                            _i--;
                            _q = qN;
                            _z += _kz * ds;
                            if (_i <= _imin) _phase = Phase::UpOutwards;
                        }
                        else
                        {
                            double ds = dsz;
                            setSegment(m, ds);
                            _k++;
                            if (_k >= _grid->_Nz)
                                setState(State::Outside);
                            else
                            {
                                _q += _kq * ds;
                                _z = zN;
                            }
                        }
                        return true;
                    }
                    case Phase::UpOutwards:
                    {
                        double RN = _grid->_Rv[_i + 1];
                        double qN = sqrt((RN - _p) * (RN + _p));
                        double zN = _grid->_zv[_k + 1];

                        int m = _grid->index(_i, 0, _k);
                        double dsq = (qN - _q) / _kq;
                        double dsz = (zN - _z) / _kz;
                        if (dsq < dsz)
                        {
                            double ds = dsq;
                            setSegment(m, ds);
                            _i++;
                            if (_i >= _grid->_NR)
                                setState(State::Outside);
                            else
                            {
                                _q = qN;
                                _z += _kz * ds;
                            }
                        }
                        else
                        {
                            double ds = dsz;
                            setSegment(m, ds);
                            _k++;
                            if (_k >= _grid->_Nz)
                                setState(State::Outside);
                            else
                            {
                                _q += _kq * ds;
                                _z = zN;
                            }
                        }
                        return true;
                    }
                    case Phase::DownInwards:
                    {
                        double RN = _grid->_Rv[_i];
                        double qN = -sqrt((RN - _p) * (RN + _p));
                        double zN = _grid->_zv[_k];

                        int m = _grid->index(_i, 0, _k);
                        double dsq = (qN - _q) / _kq;
                        double dsz = (zN - _z) / _kz;
                        if (dsq < dsz)
                        {
                            double ds = dsq;
                            setSegment(m, ds);
                            _i--;
                            _q = qN;
                            _z += _kz * ds;
                            if (_i <= _imin) _phase = Phase::DownOutwards;
                        }
                        else
                        {
                            double ds = dsz;
                            setSegment(m, ds);
                            _k--;
                            if (_k < 0)
                                setState(State::Outside);
                            else
                            {
                                _q += _kq * ds;
                                _z = zN;
                            }
                        }
                        return true;
                    }
                    case Phase::DownOutwards:
                    {
                        double RN = _grid->_Rv[_i + 1];
                        double qN = sqrt((RN - _p) * (RN + _p));
                        double zN = _grid->_zv[_k];

                        int m = _grid->index(_i, 0, _k);
                        double dsq = (qN - _q) / _kq;
                        double dsz = (zN - _z) / _kz;
                        if (dsq < dsz)
                        {
                            double ds = dsq;
                            setSegment(m, ds);
                            _i++;
                            if (_i >= _grid->_NR)
                                setState(State::Outside);
                            else
                            {
                                _q = qN;
                                _z += _kz * ds;
                            }
                        }
                        else
                        {
                            double ds = dsz;
                            setSegment(m, ds);
                            _k--;
                            if (_k < 0)
                                setState(State::Outside);
                            else
                            {
                                _q += _kq * ds;
                                _z = zN;
                            }
                        }
                        return true;
                    }
                }
            }

            case State::Outside:
            {
            }
        }
        return false;
    }
};

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> Cylinder3DSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

//////////////////////////////////////////////////////////////////////

void Cylinder3DSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // cylinders
    for (double R : _Rv) outfile->writeCircle(R);

    // meridional planes
    for (double phi : _phiv)
        outfile->writeLine(_Rv[0] * cos(phi), _Rv[0] * sin(phi), _Rv[_NR] * cos(phi), _Rv[_NR] * sin(phi));
}

//////////////////////////////////////////////////////////////////////

void Cylinder3DSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    // cylinders
    for (double R : _Rv)
    {
        outfile->writeLine(R, _zv[0], R, _zv[_Nz]);
        outfile->writeLine(-R, _zv[0], -R, _zv[_Nz]);
    }

    // horizontal planes
    for (double z : _zv)
    {
        outfile->writeLine(-_Rv[_NR], z, -_Rv[0], z);
        outfile->writeLine(_Rv[0], z, _Rv[_NR], z);
    }
}

//////////////////////////////////////////////////////////////////////

void Cylinder3DSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    write_xz(outfile);
}

//////////////////////////////////////////////////////////////////////

void Cylinder3DSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // for all cylinders
    for (double R : _Rv)
    {
        // draw the intersections with the horizontal planes
        for (double z : _zv) outfile->writeCircle(R, z);

        // draw the intersections with the meridional planes
        for (double phi : _phiv)
            outfile->writeLine(R * cos(phi), R * sin(phi), _zv[0], R * cos(phi), R * sin(phi), _zv[_Nz]);
    }

    // draw the intersections of the horizontal and meridional planes
    for (double z : _zv)
    {
        for (double phi : _phiv)
        {
            outfile->writeLine(_Rv[0] * cos(phi), _Rv[0] * sin(phi), z, _Rv[_NR] * cos(phi), _Rv[_NR] * sin(phi), z);
        }
    }
}

//////////////////////////////////////////////////////////////////////

int Cylinder3DSpatialGrid::index(int i, int j, int k) const
{
    return k + (j + i * _Nphi) * _Nz;
}

//////////////////////////////////////////////////////////////////////

bool Cylinder3DSpatialGrid::getCoords(int m, double& Rmin, double& phimin, double& zmin, double& Rmax, double& phimax,
                                      double& zmax) const
{
    if (m < 0 || m >= _Ncells) return false;

    int i = m / (_Nphi * _Nz);
    int j = (m / _Nz) % _Nphi;
    int k = m % _Nz;

    Rmin = _Rv[i];
    phimin = _phiv[j];
    zmin = _zv[k];
    Rmax = _Rv[i + 1];
    phimax = _phiv[j + 1];
    zmax = _zv[k + 1];
    return true;
}

//////////////////////////////////////////////////////////////////////

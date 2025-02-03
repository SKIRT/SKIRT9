/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Cylinder3DSpatialGrid.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Quadratic.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // small value used to check for parallel directions
    constexpr double EPS = 1e-12;
}

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

    // limit the epsilon we use for progressing the path to a value smaller than the hole and/or the first bin
    _hasHole = _Rv[0] > 0.;
    _eps = min(EPS * maxRadius(), 0.1 * (_hasHole ? min(_Rv[0], _Rv[1] - _Rv[0]) : _Rv[1]));

    // if the cylinder has no hole, create an artificial hole larger than the epsilon we use for progressing the path
    // so that the segment generator has a chance to reset the phi bin index when crossing the origin
    if (!_hasHole) _Rv[0] = 2. * _eps;

    // verify that the azimuth bins are smaller than 120 degrees, with some leeway,
    // so that the segment generator never inadvertently intersects the path with the reflected phi bin border
    constexpr double maxPhi = 0.7 * M_PI;
    for (int j = 0; j != _Nphi; ++j)
        if (_phiv[j + 1] - _phiv[j] > maxPhi) throw FATALERROR("Azimuth bin is wider than 120 deg");

    // pre-calculate sines and cosines for azimuthal bin borders; make sure that the boundary values are exact
    _sinv = sin(_phiv);
    _cosv = cos(_phiv);
    _sinv[0] = 0.;
    _cosv[0] = -1.;
    _sinv[_Nphi] = 0.;
    _cosv[_Nphi] = -1.;
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
    double _eps{0.};             // small value relative to domain size
    double _kq2{0.};             // cylindrical path direction, squared
    int _i{-1}, _j{-1}, _k{-1};  // current bin indices

public:
    MySegmentGenerator(const Cylinder3DSpatialGrid* grid) : _grid(grid), _eps(grid->_eps) {}

    // determines and sets the indices i, j and k of the cell containing the current position
    //   i is set to -1 if the position is inside Rmin and to NR if the position is outside Rmax
    //   j is clipped to the range 0..Nphi-1
    //   k is set to -1 if the position is below zmin or above zmax
    // returns true if the position is inside the grid, false if it is outside
    bool setCellIndices()
    {
        double R, phi, z;
        r().cylindrical(R, phi, z);
        _i = NR::locate(_grid->_Rv, R);
        _j = NR::locateClip(_grid->_phiv, phi);
        _k = NR::locateFail(_grid->_zv, z);
        return _i < _grid->_NR && _k >= 0;
    }

    // sets the state to outside and returns false
    bool abortPath()
    {
        setState(State::Outside);
        return false;
    }

    // returns the distance to the first intersection (or 0 if there is no intersection)
    // between the current path and the cylinder with the given bin index
    double firstIntersectionCylinder(int i)
    {
        if (abs(_kq2) < EPS) return 0.;
        double b = rx() * kx() + ry() * ky();
        double c = rx() * rx() + ry() * ry() - _grid->_Rv[i] * _grid->_Rv[i];
        return Quadratic::smallestPositiveSolution(b / _kq2, c / _kq2);
    }

    // returns the distance to the intersection (or 0 if there is no intersection)
    // between the current path and the meridional plane with the given bin index
    double intersectionMeridionalPlane(int j)
    {
        double q = kx() * _grid->_sinv[j] - ky() * _grid->_cosv[j];
        if (abs(q) < EPS) return 0.;
        return -(rx() * _grid->_sinv[j] - ry() * _grid->_cosv[j]) / q;
    }

    // returns the distance to the intersection (or 0 if there is no intersection)
    // between the current path and the horizontal plane with the given bin index
    double intersectionHorizontalPlane(int k)
    {
        if (abs(kz()) < EPS) return 0.;
        return (_grid->_zv[k] - rz()) / kz();
    }

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // initialize radial path direction
                _kq2 = kx() * kx() + ky() * ky();

                // try moving the path inside the grid; if this is impossible, return an empty path
                // keep track of cumulative length of all subsegments
                double cumds = 0.;

                // --> R direction
                double R = sqrt(rx() * rx() + ry() * ry());
                if (R > _grid->maxRadius())
                {
                    double ds = firstIntersectionCylinder(_grid->_NR);
                    if (ds <= 0.) return abortPath();
                    propagater(ds + _eps);
                    cumds += ds;
                }

                // --> z direction
                if (rz() < _grid->_zv[0])
                {
                    double ds = intersectionHorizontalPlane(0);
                    if (ds <= 0.) return abortPath();
                    propagater(ds + _eps);
                    cumds += ds;
                }
                else if (rz() > _grid->_zv[_grid->_Nz])
                {
                    double ds = intersectionHorizontalPlane(_grid->_Nz);
                    if (ds <= 0.) return abortPath();
                    propagater(ds + _eps);
                    cumds += ds;
                }

                // determine the grid cell we are in; abort in case of numerical inaccuracies
                if (!setCellIndices()) return abortPath();
                setState(State::Inside);

                // if the position was outside the cylinder, return an empty path segment with the cumulative length
                if (cumds)
                {
                    setEmptySegment(cumds);
                    return true;
                }

                // otherwise fall through to determine the first actual segment
            }

            // intentionally falls through
            case State::Inside:
            {
                // if we're not inside the real or artificial hole, proceed to the next boundary in the regular way
                if (_i >= 0)
                {
                    // remember the indices of the current cell
                    int icur = _i;
                    int jcur = _j;
                    int kcur = _k;

                    // calculate the distance travelled inside the cell by considering
                    // the potential exit points for each of the six cell boundaries;
                    // the smallest positive intersection "distance" wins.
                    double ds = DBL_MAX;  // very large, but not infinity (so that infinite values are discarded)

                    // inner radial boundary (always nonzero)
                    {
                        double s = firstIntersectionCylinder(icur);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur - 1;  // may be decremented to -1 (inside the innermost boundary)
                            _j = jcur;
                            _k = kcur;
                        }
                    }

                    // outer radial boundary
                    {
                        double s = firstIntersectionCylinder(icur + 1);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur + 1;  // may be incremented to NR (beyond the outermost boundary)
                            _j = jcur;
                            _k = kcur;
                        }
                    }

                    // clockwise azimuthal boundary
                    {
                        double s = intersectionMeridionalPlane(jcur);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur > 0 ? jcur - 1 : _grid->_Nphi - 1;  //scroll from -pi to pi
                            _k = kcur;
                        }
                    }

                    // anticlockwise azimuthal boundary
                    {
                        double s = intersectionMeridionalPlane(jcur + 1);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = (jcur + 1) % _grid->_Nphi;  //scroll from pi to -pi
                            _k = kcur;
                        }
                    }

                    // lower horizontal boundary
                    {
                        double s = intersectionHorizontalPlane(kcur);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur;
                            _k = kcur - 1;  // may be decremented to -1 (beyond the lower boundary)
                        }
                    }

                    // upper horizontal boundary
                    {
                        double s = intersectionHorizontalPlane(kcur + 1);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur;
                            _k = kcur + 1;  // may be decremented to Nz (beyond the upper boundary)
                        }
                    }

                    // if no exit point was found, abort the path
                    if (ds == DBL_MAX) return abortPath();

                    // add a segment to the path
                    setSegment(_grid->index(icur, jcur, kcur), ds);
                    propagater(ds + _eps);
                    if (_i >= _grid->_NR || _k < 0 || _k >= _grid->_Nz) setState(State::Outside);
                }

                // if we're inside the artificial hole, stay in the same azimuthal bin but allow moving
                // vertically until hitting the artificial hole radius;
                // then recalculate the bin indices (the phi bin index changes when crossing the origin)
                else if (!_grid->_hasHole)
                {
                    // remember the vertical index
                    int kcur = _k;

                    // calculate the distance travelled inside the hole by considering the hole boundary
                    // and the horizontal planes; the smallest positive intersection "distance" wins.
                    double ds = DBL_MAX;  // very large, but not infinity (so that infinite values are discarded)
                    {
                        double s = firstIntersectionCylinder(0);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                        }
                    }
                    {
                        double s = intersectionHorizontalPlane(kcur);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _k = kcur - 1;  // may be decremented to -1 (beyond the lower boundary)
                        }
                    }
                    {
                        double s = intersectionHorizontalPlane(kcur + 1);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _k = kcur + 1;  // may be decremented to Nz (beyond the upper boundary)
                        }
                    }
                    // if no exit point was found, abort the path
                    if (ds == DBL_MAX) return abortPath();

                    // add a regular segment to the path in the original azimuth bin
                    setSegment(_grid->index(0, _j, kcur), ds);
                    propagater(ds + _eps);

                    // if we exit through the hole boundary, recalculate the bin indices
                    if (_k == kcur)
                    {
                        if (!setCellIndices()) return abortPath();
                    }
                    // otherwise we may have exited through the roof or floor
                    else if (_k < 0 || _k >= _grid->_Nz)
                    {
                        setState(State::Outside);
                    }
                }

                // if we're inside the real hole, skip to the hole radius in one empty segment step
                // and recalculate the bin indices (the phi bin index changes when crossing the origin)
                else
                {
                    double ds = firstIntersectionCylinder(0);
                    if (ds <= 0.) return abortPath();
                    setEmptySegment(ds);
                    propagater(ds + _eps);
                    if (!setCellIndices()) return abortPath();
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

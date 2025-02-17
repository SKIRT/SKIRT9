/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Sphere3DSpatialGrid.hpp"
#include "Cubic.hpp"
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

void Sphere3DSpatialGrid::setupSelfAfter()
{
    SphereSpatialGrid::setupSelfAfter();

    // ---- radial ----

    // set up the radial grid
    _Nr = _meshRadial->numBins();
    _rv = _meshRadial->mesh() * (maxRadius() - minRadius()) + minRadius();

    // limit the epsilon we use for progressing the path to a value smaller than the hole and/or the first bin
    bool hasHole = _rv[0] > 0.;
    _eps = min(EPS * maxRadius(), 0.1 * (hasHole ? min(_rv[0], _rv[1] - _rv[0]) : _rv[1]));

    // if the cylinder has no hole, create an artificial hole larger than the epsilon we use for progressing the path
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

    // ---- nr of cells ----

    // cash the final number of cells
    _Ncells = _Nr * _Ntheta * _Nphi;
}

//////////////////////////////////////////////////////////////////////

int Sphere3DSpatialGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int Sphere3DSpatialGrid::numCells() const
{
    return _Ncells;
}

//////////////////////////////////////////////////////////////////////

double Sphere3DSpatialGrid::volume(int m) const
{
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(m, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        return (1. / 3.) * Cubic::pow3(rmin, rmax) * (cos(thetamin) - cos(thetamax)) * (phimax - phimin);
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double Sphere3DSpatialGrid::diagonal(int m) const
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

int Sphere3DSpatialGrid::cellIndex(Position bfr) const
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

Position Sphere3DSpatialGrid::centralPositionInCell(int m) const
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

Position Sphere3DSpatialGrid::randomPositionInCell(int m) const
{
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(m, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        double r = cbrt(Cubic::pow3(rmin) + Cubic::pow3(rmin, rmax) * random()->uniform());
        double theta = acos(cos(thetamin) + (cos(thetamax) - cos(thetamin)) * random()->uniform());
        double phi = phimin + (phimax - phimin) * random()->uniform();
        return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

class Sphere3DSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const Sphere3DSpatialGrid* _grid{nullptr};
    double _eps{0.};             // small value relative to domain size
    int _i{-1}, _j{-1}, _k{-1};  // bin indices

public:
    MySegmentGenerator(const Sphere3DSpatialGrid* grid) : _grid(grid), _eps(grid->_eps) {}

    // determines and sets the indices i, j and k of the cell containing the current position
    //   i is set to -1 if the position is inside rmin and to Nr if the position is outside rmax
    //   j is clipped to the range 0..Ntheta-1
    //   k is clipped to the range 0..Nphi-1
    // returns true if the position is inside rmax, false if it is outside rmax
    bool setCellIndices()
    {
        double radius, theta, phi;
        r().spherical(radius, theta, phi);
        _i = NR::locate(_grid->_rv, radius);
        _j = NR::locateClip(_grid->_thetav, theta);
        _k = NR::locateClip(_grid->_phiv, phi);
        return _i < _grid->_Nr;
    }

    // sets the state to outside and returns false
    bool abortPath()
    {
        setState(State::Outside);
        return false;
    }

    // returns the distance to the first intersection (or 0 if there is no intersection)
    // the current path and the sphere with given bin index
    // or 0 if there is no intersection
    double firstIntersectionSphere(int i)
    {
        return Quadratic::smallestPositiveSolution(Vec::dot(r(), k()), r().norm2() - _grid->_rv[i] * _grid->_rv[i]);
    }

    // returns the smallest positive solution of a*x^2 + 2*b*x + c = 0, or 0 if there is no positive solution
    static double smallestPositiveSolution(double a, double b, double c)
    {
        if (abs(a) > EPS) return Quadratic::smallestPositiveSolution(b / a, c / a);
        double x = -0.5 * c / b;
        if (x > 0.) return x;
        return 0.;
    }

    // returns the distance to the first intersection (or 0 if there is no intersection)
    // between the current path and the cone with given bin index
    // (the degenarate cone with zero cosine is treated separately)
    double firstIntersectionCone(int j)
    {
        double c = _grid->_cv[j];
        return c ? smallestPositiveSolution(c * c - kz() * kz(), c * c * Vec::dot(r(), k()) - rz() * kz(),
                                            c * c * r().norm2() - rz() * rz())
                 : -rz() / kz();  // degenerate cone identical to xy-plane
    }

    // returns the distance to the intersection (or 0 if there is no intersection)
    // between the current path and the meridional plane with the given bin index
    double intersectionMeridionalPlane(int k)
    {
        double q = kx() * _grid->_sinv[k] - ky() * _grid->_cosv[k];
        if (abs(q) < EPS) return 0.;
        return -(rx() * _grid->_sinv[k] - ry() * _grid->_cosv[k]) / q;
    }

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // if necessary, try moving the path inside the grid
                if (r().norm() > _grid->maxRadius())
                {
                    // find intersection; abort if there is none
                    double ds = firstIntersectionSphere(_grid->_Nr);
                    if (ds <= 0.) return abortPath();

                    // propagate the path to the intersection; abort in case of numerical inaccuracies
                    propagater(ds + _eps);
                    if (!setCellIndices()) return abortPath();

                    // return an empty path segment with the appropriate length
                    setEmptySegment(ds);
                    setState(State::Inside);
                    return true;
                }

                // the original position was inside the grid
                if (!setCellIndices()) return abortPath();  // abort in case of numerical inaccuracies
                setState(State::Inside);

                // fall through to determine the first actual segment
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
                        double s = firstIntersectionSphere(icur);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur - 1;  // may be decremented to -1 (inside the innermost boundary)
                            _j = jcur;
                            _k = kcur;
                        }
                    }

                    // outer radial boundary
                    {
                        double s = firstIntersectionSphere(icur + 1);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur + 1;  // may be incremented to Nr (beyond the outermost boundary)
                            _j = jcur;
                            _k = kcur;
                        }
                    }

                    // upper angular boundary (not applicable to uppermost cell)
                    if (jcur > 0)
                    {
                        double s = firstIntersectionCone(jcur);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur - 1;
                            _k = kcur;
                        }
                    }

                    // lower angular boundary (not applicable to lowest cell)
                    if (jcur < _grid->_Ntheta - 1)
                    {
                        double s = firstIntersectionCone(jcur + 1);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur + 1;
                            _k = kcur;
                        }
                    }

                    // clockwise azimuthal boundary
                    {
                        double s = intersectionMeridionalPlane(kcur);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur;
                            _k = kcur > 0 ? kcur - 1 : _grid->_Nphi - 1;  //scroll from -pi to pi
                        }
                    }

                    // anticlockwise azimuthal boundary
                    {
                        double s = intersectionMeridionalPlane(kcur + 1);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur;
                            _k = (kcur + 1) % _grid->_Nphi;  //scroll from pi to -pi
                        }
                    }

                    // if no exit point was found, abort the path
                    if (ds == DBL_MAX) return abortPath();

                    // add a segment to the path
                    setSegment(_grid->index(icur, jcur, kcur), ds);
                    propagater(ds + _eps);
                    if (_i >= _grid->_Nr) setState(State::Outside);
                }

                // if we're inside the hole, skip to the hole radius in one empty segment step
                // and recalculate the bin indices (the phi bin index changes when crossing the origin)
                else
                {
                    double ds = firstIntersectionSphere(0);
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

std::unique_ptr<PathSegmentGenerator> Sphere3DSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

//////////////////////////////////////////////////////////////////////

void Sphere3DSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // spheres
    for (double r : _rv) outfile->writeCircle(r);

    // meridional planes
    for (double phi : _phiv)
        outfile->writeLine(_rv[0] * cos(phi), _rv[0] * sin(phi), _rv[_Nr] * cos(phi), _rv[_Nr] * sin(phi));
}

//////////////////////////////////////////////////////////////////////

void Sphere3DSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
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

void Sphere3DSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    write_xz(outfile);
}

//////////////////////////////////////////////////////////////////////

void Sphere3DSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
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

int Sphere3DSpatialGrid::index(int i, int j, int k) const
{
    return k + (j + i * _Ntheta) * _Nphi;
}

//////////////////////////////////////////////////////////////////////

bool Sphere3DSpatialGrid::getCoords(int m, double& rmin, double& thetamin, double& phimin, double& rmax,
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

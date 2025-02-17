/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Sphere2DSpatialGrid.hpp"
#include "Cubic.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Quadratic.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void Sphere2DSpatialGrid::setupSelfAfter()
{
    SphereSpatialGrid::setupSelfAfter();

    // set up the radial grid
    _Nr = _meshRadial->numBins();
    _rv = _meshRadial->mesh() * (maxRadius() - minRadius()) + minRadius();

    // set up the polar grid; pre-calculate the cosines for each angular boundary
    _Ntheta = initPolarGrid(_thetav, _cv, _meshPolar);

    // cash the final number of cells
    _Ncells = _Nr * _Ntheta;
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::dimension() const
{
    return 2;
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::numCells() const
{
    return _Ncells;
}

//////////////////////////////////////////////////////////////////////

double Sphere2DSpatialGrid::volume(int m) const
{
    double rmin, thetamin, rmax, thetamax;
    if (getCoords(m, rmin, thetamin, rmax, thetamax))
    {
        return (2. / 3.) * M_PI * Cubic::pow3(rmin, rmax) * (cos(thetamin) - cos(thetamax));
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double Sphere2DSpatialGrid::diagonal(int m) const
{
    double rmin, thetamin, rmax, thetamax;
    if (getCoords(m, rmin, thetamin, rmax, thetamax))
    {
        Position p0(rmin, thetamin, 0., Position::CoordinateSystem::SPHERICAL);
        Position p1(rmax, thetamax, 0., Position::CoordinateSystem::SPHERICAL);
        return (p1 - p0).norm();
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::cellIndex(Position bfr) const
{
    double r, theta, phi;
    bfr.spherical(r, theta, phi);

    int i = NR::locateFail(_rv, r);
    if (i < 0) return -1;
    int j = NR::locateClip(_thetav, theta);

    return index(i, j);
}

//////////////////////////////////////////////////////////////////////

Position Sphere2DSpatialGrid::centralPositionInCell(int m) const
{
    double rmin, thetamin, rmax, thetamax;
    if (getCoords(m, rmin, thetamin, rmax, thetamax))
    {
        double r = 0.5 * (rmin + rmax);
        double theta = 0.5 * (thetamin + thetamax);
        double phi = 0.;
        return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

Position Sphere2DSpatialGrid::randomPositionInCell(int m) const
{
    double rmin, thetamin, rmax, thetamax;
    if (getCoords(m, rmin, thetamin, rmax, thetamax))
    {
        double r = cbrt(Cubic::pow3(rmin) + Cubic::pow3(rmin, rmax) * random()->uniform());
        double theta = acos(cos(thetamin) + (cos(thetamax) - cos(thetamin)) * random()->uniform());
        double phi = 2.0 * M_PI * random()->uniform();
        return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // returns the smallest positive solution of a*x^2 + 2*b*x + c = 0, or 0 if there is no positive solution
    double smallestPositiveSolution(double a, double b, double c)
    {
        if (fabs(a) > 1e-9) return Quadratic::smallestPositiveSolution(b / a, c / a);
        double x = -0.5 * c / b;
        if (x > 0.) return x;
        return 0.;
    }

    // returns the distance to the first intersection between the ray (bfr,bfk) and the sphere with given radius,
    // or 0 if there is no intersection
    double firstIntersectionSphere(Vec bfr, Vec bfk, double r)
    {
        return Quadratic::smallestPositiveSolution(Vec::dot(bfr, bfk), bfr.norm2() - r * r);
    }

    // returns the distance to the first intersection between the ray (bfr,bfk) and the cone with given cos(theta),
    // or 0 if there is no intersection (the degenarate cone with zero cosine is treated separately)
    double firstIntersectionCone(Vec bfr, Vec bfk, double c)
    {
        return c ? smallestPositiveSolution(c * c - bfk.z() * bfk.z(), c * c * Vec::dot(bfr, bfk) - bfr.z() * bfk.z(),
                                            c * c * bfr.norm2() - bfr.z() * bfr.z())
                 : -bfr.z() / bfk.z();  // degenerate cone identical to xy-plane
    }
}

//////////////////////////////////////////////////////////////////////

class Sphere2DSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const Sphere2DSpatialGrid* _grid{nullptr};
    double _eps{0.};     // small value relative to domain size
    int _i{-1}, _j{-1};  // bin indices

public:
    MySegmentGenerator(const Sphere2DSpatialGrid* grid) : _grid(grid), _eps(1e-11 * grid->maxRadius()) {}

    // determines and sets the indices i and j of the cell containing the current position
    //   i is set to -1 if the position is inside rmin and to Nr if the position is outside rmax
    //   j is clipped to the range 0..Ntheta-1
    // returns true if the position is inside rmax, false if it is outside rmax
    bool setCellIndices()
    {
        double radius, theta, phi;
        r().spherical(radius, theta, phi);
        _i = NR::locate(_grid->_rv, radius);
        _j = NR::locateClip(_grid->_thetav, theta);
        return _i < _grid->_Nr;
    }

    // sets the state to outside and returns false
    bool abortPath()
    {
        setState(State::Outside);
        return false;
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
                    double ds = firstIntersectionSphere(r(), k(), _grid->maxRadius());
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
                // if necessary, move it away from the origin so that it has meaningful cell indices
                if (r().norm() < _eps) propagater(_eps);
                if (!setCellIndices()) return abortPath();  // abort in case of numerical inaccuracies
                setState(State::Inside);

                // fall through to determine the first actual segment
            }

            // intentionally falls through
            case State::Inside:
            {
                while (true)  // the loop is executed more than once only if no exit point is found
                {
                    // remember the indices of the current cell
                    int icur = _i;
                    int jcur = _j;

                    // calculate the distance travelled inside the cell by considering
                    // the potential exit points for each of the four cell boundaries;
                    // the smallest positive intersection "distance" wins.
                    double ds = DBL_MAX;  // very large, but not infinity (so that infinite values are discarded)

                    // inner radial boundary (not applicable to innermost cell if its radius is zero)
                    if (icur > 0 || (icur == 0 && _grid->_rv[0] > 0.))
                    {
                        double s = firstIntersectionSphere(r(), k(), _grid->_rv[icur]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur - 1;  // may be decremented to -1 (inside the innermost boundary)
                            _j = jcur;
                        }
                    }

                    // outer radial boundary
                    {
                        double s = firstIntersectionSphere(r(), k(), _grid->_rv[icur + 1]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur + 1;  // may be incremented to Nr (beyond the outermost boundary)
                            _j = jcur;
                        }
                    }

                    // upper angular boundary (not applicable to uppermost cell)
                    if (jcur > 0)
                    {
                        double s = firstIntersectionCone(r(), k(), _grid->_cv[jcur]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur - 1;
                        }
                    }

                    // lower angular boundary (not applicable to lowest cell)
                    if (jcur < _grid->_Ntheta - 1)
                    {
                        double s = firstIntersectionCone(r(), k(), _grid->_cv[jcur + 1]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur + 1;
                        }
                    }

                    // if an exit point was found, add a segment to the path
                    if (_i != icur || _j != jcur)
                    {
                        setSegment(_grid->index(icur, jcur), ds);
                        propagater(ds + _eps);
                        if (_i >= _grid->_Nr) setState(State::Outside);
                        return true;
                    }

                    // otherwise, move a tiny bit along the path and reset the current cell indices
                    // if the new current point is outside the grid, there is no segment to return
                    propagater(_eps);
                    if (!setCellIndices()) return abortPath();

                    // try again from the start of this loop
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

std::unique_ptr<PathSegmentGenerator> Sphere2DSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

//////////////////////////////////////////////////////////////////////

void Sphere2DSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _Nr; i++) outfile->writeCircle(_rv[i]);
}

//////////////////////////////////////////////////////////////////////

void Sphere2DSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    // spheres
    for (double r : _rv) outfile->writeCircle(r);

    // inclined planes
    for (double theta : _thetav)
    {
        outfile->writeLine(_rv[0] * sin(theta), _rv[0] * cos(theta), _rv[_Nr] * sin(theta), _rv[_Nr] * cos(theta));
        outfile->writeLine(-_rv[0] * sin(theta), -_rv[0] * cos(theta), -_rv[_Nr] * sin(theta), -_rv[_Nr] * cos(theta));
    }
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::index(int i, int j) const
{
    return j + _Ntheta * i;
}

//////////////////////////////////////////////////////////////////////

bool Sphere2DSpatialGrid::getCoords(int m, double& rmin, double& thetamin, double& rmax, double& thetamax) const
{
    if (m < 0 || m >= _Ncells) return false;

    int i = m / _Ntheta;
    int j = m % _Ntheta;

    rmin = _rv[i];
    thetamin = _thetav[j];
    rmax = _rv[i + 1];
    thetamax = _thetav[j + 1];
    return true;
}

//////////////////////////////////////////////////////////////////////

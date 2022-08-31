/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Sphere2DSpatialGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPath.hpp"
#include "SpatialGridPlotFile.hpp"

//////////////////////////////////////////////////////////////////////

void Sphere2DSpatialGrid::setupSelfAfter()
{
    // Set up the radial grid
    _Nr = _meshRadial->numBins();
    _rv = _meshRadial->mesh() * maxRadius();

    // Set up the polar grid; pre-calculate the cosines for each angular boundary
    _Ntheta = _meshPolar->numBins();
    _thetav = _meshPolar->mesh() * M_PI;
    _cv = cos(_thetav);

    // Make sure that the boundary cosine values are exact
    _cv[0] = 1.;
    _cv[_Ntheta] = -1.;

    // The path() function in this class requires that there is a grid point corresponding to pi/2 (i.e. the xy-plane)

    // If there is a cosine value close to zero, make it exactly zero so that the test in path() succeeds
    int numzeroes = 0;
    for (int k = 1; k < _Ntheta; k++)
    {
        if (fabs(_cv[k]) < 1e-9)
        {
            numzeroes++;
            _cv[k] = 0.;
        }
    }
    if (numzeroes > 1) throw FATALERROR("There are multiple grid points very close to pi/2");

    // If there is no cosine value close to zero, add an extra grid point
    if (numzeroes == 0)
    {
        // Make temporary copy of the original grid
        Array or_thetav = _thetav;
        Array or_cv = _cv;

        // Resize the grid to contain one extra bin; this clears all values
        _Ntheta++;
        _thetav.resize(_Ntheta + 1);
        _cv.resize(_Ntheta + 1);

        // Initialize the new grid with values corresponding to the xy-plane so we can skip that index while copying
        _thetav = M_PI_2;

        // Copy the values from the original to the new grid, skipping the xy-plane
        for (int k = 0; k < _Ntheta; k++)
        {
            int target = (or_cv[k] > 0) ? k : k + 1;
            _thetav[target] = or_thetav[k];
            _cv[target] = or_cv[k];
        }
    }

    // base class setupSelfAfter() depends on initialization performed above
    SphereSpatialGrid::setupSelfAfter();
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::dimension() const
{
    return 2;
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::numCells() const
{
    return _Nr * _Ntheta;
}

//////////////////////////////////////////////////////////////////////

double Sphere2DSpatialGrid::volume(int m) const
{
    int i, k;
    invertIndex(m, i, k);
    if (i < 0 || i >= _Nr || k < 0 || k >= _Ntheta)
        return 0.0;
    else
        return (2.0 / 3.0) * M_PI * (pow(_rv[i + 1], 3) - pow(_rv[i], 3)) * (cos(_thetav[k]) - cos(_thetav[k + 1]));
}

//////////////////////////////////////////////////////////////////////

double Sphere2DSpatialGrid::diagonal(int m) const
{
    int i, k;
    invertIndex(m, i, k);
    if (i < 0 || i >= _Nr || k < 0 || k >= _Ntheta) return 0.;
    Position p1(_rv[i + 1], _thetav[k + 1], 0., Position::CoordinateSystem::SPHERICAL);
    Position p0(_rv[i], _thetav[k], 0., Position::CoordinateSystem::SPHERICAL);
    return (p1 - p0).norm();
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::cellIndex(Position bfr) const
{
    double r, theta, phi;
    bfr.spherical(r, theta, phi);
    int i = NR::locateFail(_rv, r);
    if (i < 0) return -1;
    int k = NR::locateClip(_thetav, theta);
    return index(i, k);
}

//////////////////////////////////////////////////////////////////////

Position Sphere2DSpatialGrid::centralPositionInCell(int m) const
{
    int i, k;
    invertIndex(m, i, k);
    double r = (_rv[i] + _rv[i + 1]) / 2.0;
    double theta = (_thetav[k] + _thetav[k + 1]) / 2.0;
    double phi = 0.0;
    return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
}

//////////////////////////////////////////////////////////////////////

Position Sphere2DSpatialGrid::randomPositionInCell(int m) const
{
    int i, k;
    invertIndex(m, i, k);
    double r = cbrt(_rv[i] * _rv[i] * _rv[i]
                    + (_rv[i + 1] - _rv[i]) * (_rv[i + 1] * _rv[i + 1] + _rv[i + 1] * _rv[i] + _rv[i] * _rv[i])
                          * random()->uniform());
    double theta = acos(cos(_thetav[k]) + (cos(_thetav[k + 1]) - cos(_thetav[k])) * random()->uniform());
    double phi = 2.0 * M_PI * random()->uniform();
    return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // returns the smallest positive solution of x^2 + 2*b*x + c = 0, or 0 if there is no positive solution
    double smallestPositiveSolution(double b, double c)
    {
        // x1 == -b - sqrt(b*b-c)
        // x2 == -b + sqrt(b*b-c)
        // x1*x2 == c

        if (b * b > c)  // if discriminant is negative, there are no real solutions
        {
            if (b > 0)  // x1 is always negative; x2 is positive only if c<0
            {
                if (c < 0)
                {
                    double x1 = -b - sqrt(b * b - c);
                    return c / x1;
                }
            }
            else  // x2 is always positive; x1 is positive only if c>0
            {
                double x2 = -b + sqrt(b * b - c);
                if (c > 0)
                {
                    double x1 = c / x2;
                    if (x1 < x2) return x1;
                }
                return x2;
            }
        }
        return 0;
    }

    // returns the smallest positive solution of a*x^2 + 2*b*x + c = 0, or 0 if there is no positive solution
    double smallestPositiveSolution(double a, double b, double c)
    {
        if (fabs(a) > 1e-9) return smallestPositiveSolution(b / a, c / a);
        double x = -0.5 * c / b;
        if (x > 0) return x;
        return 0;
    }

    // returns the distance to the first intersection between the ray (bfr,bfk) and the sphere with given radius,
    // or 0 if there is no intersection
    double firstIntersectionSphere(Vec bfr, Vec bfk, double r)
    {
        return smallestPositiveSolution(Vec::dot(bfr, bfk), bfr.norm2() - r * r);
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
    double _eps{0.};
    int _i{-1}, _k{-1};

public:
    MySegmentGenerator(const Sphere2DSpatialGrid* grid) : _grid(grid) {}

    // determine the indices of the cell containing the current position
    // _i will be set to -1 if the position is outside the grid
    void setCellIndices()
    {
        double radius, theta, phi;
        r().spherical(radius, theta, phi);
        _i = NR::locateFail(_grid->_rv, radius);
        _k = NR::locateClip(_grid->_thetav, theta);
    }

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // small value relative to domain size
                double rmax = _grid->maxRadius();
                _eps = 1e-11 * rmax;

                // if necessary, try moving the photon packet inside the grid
                double r2 = r().norm2();
                if (r2 > rmax * rmax)
                {
                    double ds = firstIntersectionSphere(r(), k(), rmax);
                    if (ds > 0.)
                    {
                        propagater(ds + _eps);
                        setCellIndices();
                        if (_i >= 0)
                        {
                            // return an empty path segment with the appropriate length
                            setEmptySegment(ds);
                            setState(State::Inside);
                            return true;
                        }
                    }
                    // there is no intersection with the grid; return an empty path
                    setState(State::Outside);
                    return false;
                }

                // the original position was inside the grid
                // if necessary, move the position away from the origin so that it has meaningful cell indices
                if (r2 == 0.) propagater(_eps);
                setCellIndices();
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
                    int kcur = _k;

                    // calculate the distance travelled inside the cell by considering
                    // the potential exit points for each of the four cell boundaries;
                    // the smallest positive intersection "distance" wins.
                    double ds = DBL_MAX;  // very large, but not infinity (so that infinite values are discarded)

                    // inner radial boundary (not applicable to innermost cell)
                    if (icur > 0)
                    {
                        double s = firstIntersectionSphere(r(), k(), _grid->_rv[icur]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur - 1;
                            _k = kcur;
                        }
                    }

                    // outer radial boundary (always applicable)
                    {
                        double s = firstIntersectionSphere(r(), k(), _grid->_rv[icur + 1]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur + 1;  // may be incremented beyond the outermost boundary
                            _k = kcur;
                        }
                    }

                    // upper angular boundary (not applicable to uppermost cell)
                    if (kcur > 0)
                    {
                        double s = firstIntersectionCone(r(), k(), _grid->_cv[kcur]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _k = kcur - 1;
                        }
                    }

                    // lower angular boundary (not applicable to lowest cell)
                    if (kcur < _grid->_Ntheta - 1)
                    {
                        double s = firstIntersectionCone(r(), k(), _grid->_cv[kcur + 1]);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _k = kcur + 1;
                        }
                    }

                    // if an exit point was found, add a segment to the path
                    if (_i != icur || _k != kcur)
                    {
                        setSegment(_grid->index(icur, kcur), ds);
                        propagater(ds + _eps);
                        if (_i >= _grid->_Nr) setState(State::Outside);
                        return true;
                    }

                    // otherwise, move a tiny bit along the path and reset the current cell indices
                    propagater(_eps);
                    setCellIndices();
                    if (_i < 0)
                    {
                        // the new current point is outside the grid; there is no segment to return
                        setState(State::Outside);
                        return false;
                    }
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

/*
{
    // Small value relative to domain size
    double rmax = maxRadius();
    const double eps = 1e-11 * rmax;

    // Initialize the path
    path->clear();
    Position bfr = path->position();
    Direction bfk = path->direction();

    // Move the photon packet to the first grid cell that it will pass.
    // If it does not pass any grid cell, return an empty path.
    // Otherwise calculate the distance covered and add a segment to the path.
    double r2 = bfr.norm2();
    if (r2 > rmax * rmax)
    {
        double ds = firstIntersectionSphere(bfr, bfk, rmax);
        if (!ds) return path->clear();
        path->addSegment(-1, ds);
        bfr += bfk * (ds + eps);
    }
    // Move the position a bit away from the origin so that it has a meaningful cell number
    else if (r2 == 0)
    {
        bfr += bfk * eps;
    }

    // Determine the indices of the cell containing the starting point.
    double r, theta, phi;
    bfr.spherical(r, theta, phi);
    int i = NR::locateFail(_rv, r);
    int k = NR::locateClip(_thetav, theta);

    // Start the loop over cells/path segments until we leave the grid
    int inext = i;
    int knext = k;
    while (i < _Nr && i >= 0)
    {
        // Calculate the distance travelled inside the cell by considering
        // the potential exit points for each of the four cell boundaries;
        // the smallest positive intersection "distance" wins.
        double ds = DBL_MAX;  // very large, but not infinity (so that infinite values are discarded)

        // inner radial boundary (not applicable to innermost cell)
        if (i > 0)
        {
            double s = firstIntersectionSphere(bfr, bfk, _rv[i]);
            if (s > 0 && s < ds)
            {
                ds = s;
                inext = i - 1;
                knext = k;
            }
        }

        // outer radial boundary (always applicable)
        {
            double s = firstIntersectionSphere(bfr, bfk, _rv[i + 1]);
            if (s > 0 && s < ds)
            {
                ds = s;
                inext = i + 1;  // this will cause the loop to terminate if incremented beyond the outermost boundary
                knext = k;
            }
        }

        // upper angular boundary (not applicable to uppermost cell)
        if (k > 0)
        {
            double s = firstIntersectionCone(bfr, bfk, _cv[k]);
            if (s > 0 && s < ds)
            {
                ds = s;
                inext = i;
                knext = k - 1;
            }
        }

        // lower angular boundary (not applicable to lowest cell)
        if (k < _Ntheta - 1)
        {
            double s = firstIntersectionCone(bfr, bfk, _cv[k + 1]);
            if (s > 0 && s < ds)
            {
                ds = s;
                inext = i;
                knext = k + 1;
            }
        }

        // If an exit point was found, add a segment to the path,
        // move to the next current point, and update the cell indices
        if (inext != i || knext != k)
        {
            path->addSegment(index(i, k), ds);
            bfr += bfk * (ds + eps);
            i = inext;
            k = knext;
        }
        // Otherwise, move a tiny bit along the path and reset the current cell indices
        else
        {
            find<Log>()->warning("No exit point found from Spatial grid cell");
            bfr += bfk * eps;
            double r, theta, phi;
            bfr.spherical(r, theta, phi);
            i = NR::locateFail(_rv, r);
            k = NR::locateClip(_thetav, theta);
        }
    }
}
*/
//////////////////////////////////////////////////////////////////////

void Sphere2DSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    for (int i = 0; i <= _Nr; i++) outfile->writeCircle(_rv[i]);
}

//////////////////////////////////////////////////////////////////////

void Sphere2DSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    double rmax = maxRadius();
    for (int i = 0; i <= _Nr; i++)
    {
        outfile->writeCircle(_rv[i]);
    }
    for (int k = 0; k <= _Ntheta; k++)
    {
        double x = rmax * sin(_thetav[k]);
        double z = rmax * cos(_thetav[k]);
        outfile->writeLine(-x, -z, x, z);
    }
}

//////////////////////////////////////////////////////////////////////

int Sphere2DSpatialGrid::index(int i, int k) const
{
    return k + _Ntheta * i;
}

//////////////////////////////////////////////////////////////////////

void Sphere2DSpatialGrid::invertIndex(int m, int& i, int& k) const
{
    i = m / _Ntheta;
    k = m % _Ntheta;
}

//////////////////////////////////////////////////////////////////////

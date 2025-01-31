/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Sphere3DSpatialGrid.hpp"
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
    _Ntheta = _meshPolar->numBins();
    _thetav = _meshPolar->mesh() * M_PI;
    _cv = cos(_thetav);

    // make sure that the boundary cosine values are exact
    _cv[0] = 1.;
    _cv[_Ntheta] = -1.;

    // the path segment generator requires that there is a grid point corresponding to pi/2 (i.e. the xy-plane)

    // if there is a cosine value close to zero, make it exactly zero
    int numZeroes = 0;
    for (int j = 1; j < _Ntheta; j++)
    {
        if (fabs(_cv[j]) < 1e-9)
        {
            numZeroes++;
            _cv[j] = 0.;
        }
    }
    if (numZeroes > 1) throw FATALERROR("There are multiple grid points very close to pi/2");

    // if there is no cosine value close to zero, add an extra grid point
    if (numZeroes == 0)
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
        for (int j = 0; j < _Ntheta; j++)
        {
            int target = (or_cv[j] > 0) ? j : j + 1;
            _thetav[target] = or_thetav[j];
            _cv[target] = or_cv[j];
        }
    }

    // ---- azimuthal ----

    // set up the azimuthal grid
    _Nphi = _meshAzimuthal->numBins();
    _phiv = _meshAzimuthal->mesh() * (2. * M_PI) - M_PI;

    // verify that the azimuth bins are smaller than 120 degrees, with some leeway,
    // so that the segment generator never inadvertently intersects the path with the reflected phi bin border
    constexpr double maxPhi = 0.7 * M_PI;
    for (int j = 0; j != _Nphi; ++j)
        if (_phiv[j + 1] - _phiv[j] > maxPhi) throw FATALERROR("Azimuth bin is wider than 120 deg");

    // pre-calculate sines and cosines for azimuthal bin borders
    _sinv = sin(_phiv);
    _cosv = cos(_phiv);

    // make sure that the boundary values are exact
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
        return (1. / 3.) * pow3(rmin, rmax) * (cos(thetamin) - cos(thetamax)) * (phimax - phimin);
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
        double r = cbrt(pow3(rmin) + pow3(rmin, rmax) * random()->uniform());
        double theta = acos(cos(thetamin) + (cos(thetamax) - cos(thetamin)) * random()->uniform());
        double phi = phimin + (phimax - phimin) * random()->uniform();
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

class Sphere3DSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const Sphere3DSpatialGrid* _grid{nullptr};
    double _eps{0.};     // small value relative to domain size
    int _i{-1}, _j{-1};  // bin indices

public:
    MySegmentGenerator(const Sphere3DSpatialGrid* grid) : _grid(grid), _eps(1e-11 * grid->maxRadius()) {}

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
                        setSegment(_grid->index(icur, 0, jcur), ds);
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

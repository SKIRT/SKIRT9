/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Sphere3DSpatialGrid.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Quadrics.hpp"

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

    // returns the distance to the first intersection (or 0 if there is no intersection) between
    // the current path and the sphere with given bin index
    double firstIntersectionSphere(int i) { return Quadrics::firstIntersectionSphere(r(), k(), _grid->_rv[i]); }

    // returns the distance to the first intersection (or 0 if there is no intersection) between
    // the current path and the cone with given bin index
    double firstIntersectionCone(int j) { return Quadrics::firstIntersectionCone(r(), k(), _grid->_cv[j]); }

    // returns the distance to the intersection (or 0 if there is no intersection) between the
    // current path and the meridional plane with the given bin index
    double intersectionMeridionalPlane(int kbin)
    {
        return Quadrics::intersectionMeridionalPlane(r(), k(), _grid->_sinv[kbin], _grid->_cosv[kbin]);
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

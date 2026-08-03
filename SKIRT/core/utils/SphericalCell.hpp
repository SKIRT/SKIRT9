/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERICALCELL_HPP
#define SPHERICALCELL_HPP

#include "Vec.hpp"
class Box;
class Position;

//////////////////////////////////////////////////////////////////////

/** SphericalCell is a low-level class for working with basic three-dimensional cells in spherical
    coordinates. Each SphericalCell instance represents a cell bordered by:

    - two spheres centered on the origin defined by radii \f$0 \le r_\text{min} \le
    r_\text{max}\f$,

    - two polar half-cones defined by inclination angles \f$0 \le \theta_\text{min} \le
    \theta_\text{max} \le \pi\f$,

    - two meridional half-planes (with \f$r>0\f$) defined by azimuth angles \f$-\pi \le
    \varphi_\text{min} \le \varphi_\text{max} \le \pi\f$ with \f$\varphi_\text{max} -
    \varphi_\text{min} \le \pi\f$,

    \note Because of the limitations on the range of \f$\varphi\f$, a SphericalCell cannot straddle
    the negative x-axis of the Cartesian model coordinate system, and it cannot span more than half
    of the azimuth circle. Also, because of the limitations on the range of \f$\theta\f$, a
    SphericalCell cannot straddle the z-axis of the Cartesian model coordinate system.

    The class offers functions to retrieve various basic properties of the cell, such as its border
    coordinates and its volume, and for geometric operations such as determining whether a given
    Cartesian position is inside the cell or calculating the intersection with a ray. */
class SphericalCell
{
public:
    /** The default constructor creates an empty cell at the origin, i.e. it initializes all border
        coordinates to zero. */
    SphericalCell() {}

    /** This constructor initializes the cell border coordinates to the values provided as
        arguments. It does not verify that these values conform to the limits described in the
        class header. Non-conforming values lead to undefined behavior. */
    SphericalCell(double rmin, double thetamin, double phimin, double rmax, double thetamax, double phimax);

    /** This function stores the cell border coordinates in the provided arguments. */
    void extent(double& rmin, double& thetamin, double& phimin, double& rmax, double& thetamax, double& phimax) const
    {
        rmin = _rmin;
        thetamin = _thetamin;
        phimin = _phimin;
        rmax = _rmax;
        thetamax = _thetamax;
        phimax = _phimax;
    }

    /** This function returns the \f$r_\text{min}\f$ border coordinate of the cell. */
    double rmin() const { return _rmin; }

    /** This function returns the \f$\theta_\text{min}\f$ border coordinate of the cell. */
    double thetamin() const { return _thetamin; }

    /** This function returns the \f$\varphi_\text{min}\f$ border coordinate of the cell. */
    double phimin() const { return _phimin; }

    /** This function returns the \f$R_\text{max}\f$ border coordinate of the cell. */
    double rmax() const { return _rmax; }

    /** This function returns the \f$\theta_\text{max}\f$ border coordinate of the cell. */
    double thetamax() const { return _thetamax; }

    /** This function returns the \f$\varphi_\text{max}\f$ border coordinate of the cell. */
    double phimax() const { return _phimax; }

    /** This function returns the volume of the cell, given by \f$\frac{1}{3}
        \left(r_\text{max}^3-r_\text{min}^3\right)
        \left(\cos\theta_\text{min}-\cos\theta_\text{max}\right)
        \left(\varphi_\text{max}-\varphi_\text{min}\right)\f$. */
    double volume() const;

    /** This function returns the "center" of the cell in Cartesian coordinates. This position is
        defined as the halfway point between the cell borders in spherical coordinates, i.e.
        \f[\begin{aligned} r_\text{ctr} = (r_\text{min} + r_\text{max})/2 \\ \theta_\text{ctr} =
        (\theta_\text{min} + \theta_\text{max})/2 \\ \varphi_\text{ctr} = (\varphi_\text{min} +
        \varphi_\text{max})/2. \end{aligned}\f]

        As stated above the function returns this position after converting it to Cartesian
        coordinates. */
    Position center() const;

    /** This function returns true if the Cartesian position \f$\bf{r}=(x,y,z)\f$ is inside the
        cell, and false otherwise. A position on an edge or face on the "lower" side of the cell is
        considered to be contained in the cell, while a position on an edge or face on the "upper"
        side of the cell is considered \em not to be contained in the cell. This approach avoids
        duplicate containment of adjacent cells. */
    bool contains(Vec bfr) const;

    /** This function returns the Cartesian bounding box of the cell, in other words the smallest
        cuboid lined up with the Cartesian coordinate axes that encloses the cell.

        The intersection of a conical cell boundary and a spherical cell boundary is a horizontal
        circle segment. This means that the z coordinate of the intersection does not depend on
        \f$\varphi\f$, and the vertical bounds of the cell can be determined by considering the
        cell's corner points in the meridional plane.

        For the projection on the equatorial plane, things are more complicated. If the cell
        straddles an x or y axis in the azimuthal direction, or if it straddles the xy plane in the
        polar direction, the projection of the outer sphere boundary on the x and/or y axes extends
        beyond the projection of the cell corner points. The bounding box must thus enclose the
        corresponding extreme points in addition to the cell corner points. */
    Box boundingBox() const;

    /** This function intersects the receiving cell with a ray (half-line) defined by the specified
        starting position \f$\bf{r}\f$ and direction \f$\bf{k}\f$. If the ray does not intersect
        the cell, the function returns zero. If the ray does intersect the cell, the function
        returns the length of the intersection segment, or if applicable, the sum of the two
        intersection segments (because a spherical cell is concave, a ray can intersect it multiple
        times). If the starting position is inside the cell, only the portion of the ray inside the
        cell is taken into account.

        A ray that touches the cell border in a single point is not considered to intersect. A ray
        along an edge or face on the "lower" side of the cell is considered to intersect, while ray
        along an edge or face on the "upper" side of the cell is considered \em not to intersect.
        This approach avoids duplicate intersection of adjacent cells.

        <em>Implementation</em>

        The function uses the following strategy. First it calculates all intersection points with
        the bordering surfaces and sorts them on distance along the ray. Then, for each segment
        between consecutive intersections beyond the ray's starting point, it determines whether
        the segment's midpoint is inside the cell or not, and accumulates the lengths of the inside
        segments. This automatically takes care of the cases where the ray lies in one of the
        bordering surfaces.

        The ray can be represented by its parameter equation \f${\bf{x}}={\bf{r}}+s\,{\bf{k}}\f$,
        with \f${\bf{k}}\f$ a unit vector.

        The two intersection points with a radial boundary sphere \f${\bf{x}}^2=r_*^2\f$ are
        obtained by solving the quadratic equation \f$s^2 + 2\,({\bf{r}}\cdot{\bf{k}})\,s +
        ({\bf{r}}^2-r_*^2)=0\f$ for \f$s\f$.

        The two intersection points with an angular boundary cone \f$x_z^2=c^2\,{\bf{x}}^2\f$ (with
        \f$c=\cos\theta_*\f$) are obtained by solving the quadratic equation \f$(c^2-k_z^2)\,s^2 +
        2\,(c^2\,{\bf{r}}\cdot{\bf{k}}-r_z k_z)\,s + (c^2\,{\bf{r}}^2-r_z^2)=0\f$ for \f$s\f$.

        The intersection point with a meriodonal plane \f$\sin\varphi_* x = \cos\varphi_* y\f$ is
        obtained by \f[ s = -\;\frac{r_\text{x}\sin\varphi_* - r_\text{y}\cos\varphi_*}
        {k_\text{x}\sin\varphi_* - k_\text{y}\cos\varphi_*} \f] */
    double intersection(Vec r, const Vec k) const;

private:
    // These data members represent the spherical border coordinates
    double _rmin{0}, _thetamin{0}, _phimin{0}, _rmax{0}, _thetamax{0}, _phimax{0};

    // These data members remember frequently used and expensive to calculate values
    double _cosphimin{1}, _sinphimin{0}, _cosphimax{1}, _sinphimax{0};
    double _costhetamin{1}, _sinthetamin{0}, _costhetamax{1}, _sinthetamax{0};
};

//////////////////////////////////////////////////////////////////////

#endif

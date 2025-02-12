/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLCELL_HPP
#define CYLCELL_HPP

#include "Box.hpp"

//////////////////////////////////////////////////////////////////////

/** CylCell is a low-level class for working with basic three-dimensional cells in cylindrical
    coordinates. Each CylCell instance represents a cell bordered by:

    - two vertical cylinders centered on the origin defined by radii \f$0 \le R_\text{min} \le
    R_\text{max}\f$,

    - two meridional half-planes (with \f$R>0\f$) defined by azimuth angles \f$-\pi \le
    \varphi_\text{min} \le \varphi_\text{max} \le \pi\f$ with \f$\varphi_\text{max} -
    \varphi_\text{min} \le \pi\f$,

    - two horizontal planes defined by \f$z_\text{min} \le z_\text{max}\f$.

    \note Because of the limitations on the range of \f$\varphi\f$, a Cylcell cannot straddle the
    negative x-axis of the Cartesian model coordinate system, and it cannot span more than half of
    the azimuth circle.

    The class offers functions to retrieve various basic properties of the cell, such as its border
    coordinates and its volume, and for geometric operations such as determining whether a given
    Cartesian position is inside the cell. */
class CylCell
{
public:
    /** The default constructor creates an empty cell at the origin, i.e. it initializes all border
        coordinates to zero. */
    CylCell() {}

    /** This constructor initializes the cell border coordinates to the values provided as
        arguments. It does not verify that these values conform to the limits described in the
        class header. Non-comforming values lead to undefined behavior. */
    CylCell(double Rmin, double phimin, double zmin, double Rmax, double phimax, double zmax);

    /** This function stores the cell border coordinates in the provided arguments. */
    void extent(double& Rmin, double& phimin, double& zmin, double& Rmax, double& phimax, double& zmax) const
    {
        Rmin = _Rmin;
        phimin = _phimin;
        zmin = _zmin;
        Rmax = _Rmax;
        phimax = _phimax;
        zmax = _zmax;
    }

    /** This function returns the \f$R_\text{min}\f$ border coordinate of the cell. */
    double Rmin() const { return _Rmin; }

    /** This function returns the \f$\varphi_\text{min}\f$ border coordinate of the cell. */
    double phimin() const { return _phimin; }

    /** This function returns the \f$z_\text{min}\f$ border coordinate of the cell. */
    double zmin() const { return _zmin; }

    /** This function returns the \f$R_\text{max}\f$ border coordinate of the cell. */
    double Rmax() const { return _Rmax; }

    /** This function returns the \f$\varphi_\text{max}\f$ border coordinate of the cell. */
    double phimax() const { return _phimax; }

    /** This function returns the \f$z_\text{max}\f$ border coordinate of the cell. */
    double zmax() const { return _zmax; }

    /** This function returns the width \f$R_\text{max}-R_\text{min}\f$ of the cell. */
    double Rwidth() const { return _Rmax - _Rmin; }

    /** This function returns the width \f$\varphi_\text{max}-\varphi_\text{min}\f$ of the cell. */
    double phiwidth() const { return _phimax - _phimin; }

    /** This function returns the width \f$z_\text{max}-z_\text{min}\f$ of the cell. */
    double zwidth() const { return _zmax - _zmin; }

    /** This function returns the volume of the cell, given by \f$\frac{1}{2}
        (R_\text{max}^2-R_\text{min}^2) (\varphi_\text{max}-\varphi_\text{min})
        (z_\text{max}-z_\text{min})\f$. */
    double volume() const;

    /** This function returns the "center" of the cell in Cartesian coordinates. This position is
        defined as the halfway point between the cell borders in cylindrical coordinates, i.e.
        \f[\begin{aligned} R_\text{ctr} = (R_\text{min} + R_\text{max})/2 \\ \varphi_\text{ctr} =
        (\varphi_\text{min} + \varphi_\text{max})/2 \\ z_\text{ctr} = (z_\text{min} +
        z_\text{max})/2. \end{aligned}\f]

        As stated above the function returns this position after converting it to Cartesian
        coordinates. */
    Vec center() const;

    /** This function returns true if the Cartesian position \f${\bf{r}}=(x,y,z)\f$ is inside the
        cell, and false otherwise. A position on an edge or face on the "lower" side of the cell is
        considered to be contained in the cell, while a position on an edge or face on the "upper"
        side of the cell is considered \em not to be contained in the cell. This approach avoids
        duplicate containment of adjacent cells. */
    bool contains(Vec r) const;

    /** This function returns the Cartesian bounding box of the cell, in other words the smallest
        cuboid lined up with the Cartesian coordinate axes that encloses the cell.

        The bounds along the z-axis are the same for Cylindrical and Cartesian coordinates, but the
        function needs to to specifically determine the bounding rectangle projected on the xy
        plane. This bounding rectangle must of course enclose the four corner points of the cell.
        In addition, if the cell straddles one of the coordinate axes, the bounding rectangle must
        also enclose a point on that axis at radius Rmax. */
    Box boundingBox() const;

    /** This function intersects the receiving cell with a ray (half-line) defined by the specified
        starting position \f$\bf{r}\f$ and direction \f$\bf{k}\f$. If the ray does not intersect
        the cell, the function returns zero. If the ray does intersect the cell, the function
        returns the length of the intersection segment, or if applicable, the sum of the two
        intersection segments (because a cell is concave at the inner radial cylinder, a ray can
        intersect it twice). If the starting position is inside the cell, only the portion of the
        ray inside the cell is taken into account.

        A ray that touches the cell border in a single point is not considered to intersect. A ray
        along an edge or face on the "lower" side of the cell is considered to intersect, while ray
        along an edge or face on the "upper" side of the cell is considered \em not to intersect.
        This approach avoids duplicate intersection of adjacent cells.

        <em>Implementation</em>

        The function uses the following strategy. First it calculates all intersection points with
        the bordering planes and cylinders and sorts them on distance along the ray. Then, for each
        segment between consecutive intersections beyond the ray's starting point, it determines
        whether the segment's midpoint is inside the cell or not, and accumulates the lengths of
        the inside segments. This automatically takes care of the cases where the ray lies in one
        of the bordering planes or cylinders.

        The ray equation can be written as \f[\begin{cases} x = r_\text{x} + k_\text{x}s \\ y =
        r_\text{y} + k_\text{y}s \\ z = r_\text{z} + k_\text{z}s \\ \end{cases} \quad \text{with}
        \;s>0.\f]

        Intersection with a horizontal plane with equation \f$z=z_*\f$ easily yields \f$s =
        (z_*-r_\text{z})/k_\text{z}\f$.

        Intersection with a meridional plane with equation \f$\sin\varphi_* x = \cos\varphi_* y\f$
        yields \f[ s = -\;\frac{r_\text{x}\sin\varphi_* - r_\text{y}\cos\varphi_*}
        {k_\text{x}\sin\varphi_* - k_\text{y}\cos\varphi_*} \f]

        Intersection with a vertical cylinder centered on the origin with equation \f$x^2 +y^2 =
        R_*^2\f$ yields a quadratic equation of the form \f$s^2+2bs+c=0\f$ with \f[\begin{aligned}
        b &= \frac{r_\text{x}k_\text{x} + r_\text{y}k_\text{y}} {k_\text{x}^2 + k_\text{y}^2} \\ c
        &= \frac{r_\text{x}^2 + r_\text{y}^2 - R_*^2} {k_\text{x}^2 + k_\text{y}^2} \\
        \end{aligned}\f] */
    double intersection(Vec r, const Vec k) const;

private:
    // These data members represent the cylindrical border coordinates
    double _Rmin{0}, _phimin{0}, _zmin{0}, _Rmax{0}, _phimax{0}, _zmax{0};

    // These data members remember frequently used and expensive to calculate values
    double _cosphimin{1}, _sinphimin{0}, _cosphimax{1}, _sinphimax{0};
};

//////////////////////////////////////////////////////////////////////

#endif

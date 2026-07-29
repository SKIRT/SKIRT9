/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERE3DSPATIALGRID_HPP
#define SPHERE3DSPATIALGRID_HPP

#include "StructuredSphereSpatialGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** The Sphere3DSpatialGrid class is a trivial concrete subclass of the StructuredSphereSpatialGrid
    class, adding nothing beyond the structured grid itself -- see the StructuredSphereSpatialGrid
    class for a description of the grid's geometry. This class provides the path segment generator
    appropriate for a plain structured grid (i.e. one without any additional embedded geometry). */
class Sphere3DSpatialGrid : public StructuredSphereSpatialGrid
{
    ITEM_CONCRETE(Sphere3DSpatialGrid, StructuredSphereSpatialGrid, "a 3D spatial grid in spherical coordinates")
        ATTRIBUTE_TYPE_DISPLAYED_IF(Sphere3DSpatialGrid, "Level2")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a 3D spherical grid, implemented as a
        private PathSegmentGenerator subclass.

        We represent the path by its parameter equation \f${\bf{x}}={\bf{r}}+s\,{\bf{k}}\f$, and we
        assume that \f${\bf{k}}\f$ is a unit vector. The two intersection points with a radial
        boundary sphere \f${\bf{x}}^2=r_*^2\f$ are obtained by solving the quadratic equation
        \f$s^2 + 2\,({\bf{r}}\cdot{\bf{k}})\,s + ({\bf{r}}^2-r_*^2)=0\f$ for \f$s\f$.

        The two intersection points with an angular boundary cone \f$x_z^2=c^2\,{\bf{x}}^2\f$ (with
        \f$c=\cos\theta_*\f$) are obtained by solving the quadratic equation \f$(c^2-k_z^2)\,s^2 +
        2\,(c^2\,{\bf{r}}\cdot{\bf{k}}-r_z k_z)\,s + (c^2\,{\bf{r}}^2-r_z^2)=0\f$ for \f$s\f$. The
        intersection points with the reflected cone are always more distant than the other cell
        boundaries (the requirement to include the xy-plane \f$\theta=\pi/2\f$ in the grid ensures
        that this is true) and thus these phantom points are automatically ignored.

        The intersection point with a meriodonal plane \f$\sin\varphi_* x = \cos\varphi_* y\f$ is
        obtained by \f[ s = -\;\frac{r_\text{x}\sin\varphi_* - r_\text{y}\cos\varphi_*}
        {k_\text{x}\sin\varphi_* - k_\text{y}\cos\varphi_*} \f] The requirement that
        \f$\varphi_{j+1}-\varphi_j\le2\pi/3\f$ ensures that the intersection point with the
        reflected half-plane is always more distant than the other cell boundaries thus that
        phantom point is automatically ignored.

        The segment generator progresses the starting point of the path through the grid along the
        path's direction. For each step along the way, it calculates the distances to the
        intersections with all candidate borders of the current cell, and then selects the nearest
        intersection point. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

private:
    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif

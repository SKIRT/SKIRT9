/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERE1DSPATIALGRID_HPP
#define SPHERE1DSPATIALGRID_HPP

#include "Array.hpp"
#include "Mesh.hpp"
#include "SphereSpatialGrid.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** The Sphere1DSpatialGrid class is a subclass of the SphereSpatialGrid class, and represents a
    one-dimensional, spherically symmetric spatial grid. Each cell in such a grid is a spherical
    shell. A grid with \f$N_r\f$ cells is specified through a set of \f$N_r+1\f$ radial grid points
    \f$r_i, \,i=0,\ldots,N_r\f$, with \f$0\le r_\text{min} = r_0\f$, \f$r_i<r_{i+1}\f$, and
    \f$r_{N_r} = r_\text{max}\f$. */
class Sphere1DSpatialGrid : public SphereSpatialGrid
{
    ITEM_CONCRETE(Sphere1DSpatialGrid, SphereSpatialGrid, "a 1D spherically symmetric spatial grid")
        ATTRIBUTE_TYPE_ALLOWED_IF(Sphere1DSpatialGrid, "!Dimension2&!Dimension3")

        PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up a number of data members that depend on the Mesh object configured
        for this grid. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 1. */
    int dimension() const override;

    /** This function returns the number of cells \f$N_r\f$ in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. The cell index \f$m\f$
        corresponds to the radial bin with lower border index \f$i=m\f$, and the volume is easily
        calculated as \f[V = \frac{4\pi}{3}\, (r_{i+1}^3-r_i^3),\f] with \f$r_i\f$ and
        \f$r_{i+1}\f$ the inner and outer radius of the shell respectively. */
    double volume(int m) const override;

    /** This function returns the "diagonal" of the cell with index \f$m\f$. The cell index \f$m\f$
        corresponds to the radial bin with lower border index \f$i=m\f$. The function simply
        returns the cell width, i.e. the difference between the outer and inner cell radius. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. It just determines the index of the radial bin containing the position, and
        returns that number. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. The cell index
        \f$m\f$ corresponds to the radial bin with lower border index \f$i=m\f$, and the central
        radius is determined using \f[ r = \frac{r_i + r_{i+1}}{2}. \f] The returned position is
        arbitrarily located on the x-axis. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location in the cell with index \f$m\f$. The cell index
        \f$m\f$ corresponds to the radial bin with lower border index \f$i=m\f$, and a random
        radius is determined using \f[ r = \left( r_i^3 + {\cal{X}} \,(r_{i+1}^3-r_i^3)
        \right)^{1/3} \f] with \f${\cal{X}}\f$ a random deviate. This random radius is combined
        with a random position on the unit sphere to obtain a random position in the cell. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a spherical grid, implemented as a private
        PathSegmentGenerator subclass.

        \image html Sphere1DSpatialGrid.png

        The figure above is drawn in the plane containing the coordinate system origin \f$\bf{o}\f$
        and the path under consideration, defined by its starting point \f$\bf{r}\f$ and direction
        vector \f$\bf{k}\f$. Now consider the point of closest approach \f$\bf{c}\f$. The distance
        \f$q\f$ from \f$\bf{c}\f$ to \f$\bf{r}\f$ can be obtained by projecting the position vector
        \f${\bf{r}}\f$ onto the line formed by the path. Assuming \f$\bf{k}\f$ is normalized to
        unity, this distance is given by the dot product \f$q=\bf{r}\cdot\bf{k}\f$. The resulting
        value is negative if \f$\bf{r}\f$ is before \f$\bf{c}\f$ (i.e. it is going inward) and
        positive if \f$\bf{r}\f$ is after \f$\bf{c}\f$ (i.e. it is going outward). We further
        define the impact parameter \f$p\f$ as the distance of closest approach. From the
        rectangular triangles illustrated in the figure it is easily seen that \f$p^2 + q^2 =
        r^2\f$ and \f$p^2 + q_*^2 = r_*^2\f$. The segment generator uses these relations to step
        the path through the consecutive spherical bins, updating the value of \f$q\f$ along the
        way as a proxy for updating \f${\bf{r}}\f$. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    //======================== Data Members ========================

private:
    int _Nr{0};
    Array _rv;

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

////////////////////////////////////////////////////////////////////

#endif

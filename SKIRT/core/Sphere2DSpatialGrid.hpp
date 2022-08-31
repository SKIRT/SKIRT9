/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERE2DSPATIALGRID_HPP
#define SPHERE2DSPATIALGRID_HPP

#include "Array.hpp"
#include "Mesh.hpp"
#include "SphereSpatialGrid.hpp"
class Random;

//////////////////////////////////////////////////////////////////////

/** The Sphere2DSpatialGrid class is subclass of the SphereSpatialGrid class, and represents
    axisymmetric spatial grids based on spherical coordinates. The grid is defined in the
    meridional plane and rotated around the Z-axis. The meridional grid is specified through a set
    of \f$N_r+1\f$ radial grid points \f$r_i\f$ (with \f$i=0,\ldots,N_r\f$) and a set of
    \f$N_\theta+1\f$ angular grid points \f$\theta_k\f$ (with \f$k=0,\ldots,N_\theta\f$). In total
    there are \f$N_{\text{cells}} = N_r\,N_\theta\f$ cells in the grid. */
class Sphere2DSpatialGrid : public SphereSpatialGrid
{
    ITEM_CONCRETE(Sphere2DSpatialGrid, SphereSpatialGrid, "an axisymmetric spatial grid in spherical coordinates")
        ATTRIBUTE_TYPE_ALLOWED_IF(Sphere2DSpatialGrid, "!Dimension3")
        ATTRIBUTE_TYPE_DISPLAYED_IF(Sphere2DSpatialGrid, "Level2")

        PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

        PROPERTY_ITEM(meshPolar, Mesh, "the bin distribution in the polar direction")
        ATTRIBUTE_DEFAULT_VALUE(meshPolar, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up a number of data members. It also pre-calculates and stores the
        opening angle cosines \f$c_k = \cos\theta_k\f$ for the boundary cones in the angular grid,
        to help speed up calculations in the path() method. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 2 for this class. */
    int dimension() const override;

    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. In this class, the
        function determines the radial and angular bin indices \f$i\f$ and \f$k\f$ that correspond
        to the cell index \f$m\f$, and then calculates the volume as \f[ V = \frac{2\pi}{3}
        \left(r_{i+1}^3-r_i^3\right) \left(\cos\theta_k-\cos\theta_{k+1}\right). \f] */
    double volume(int m) const override;

    /** This function returns the "diagonal" of the cell with index \f$m\f$. For 2D spherical
        grids, it returns the distance between the outer/upper and inner/lower corners of the cell
        in the meridional plane. */
    double diagonal(int m) const override;

    /** This function returns the number of the cell that contains the position \f${\bf{r}}\f$. In
        this class, the function determines the radial and angular bin indices and calculates the
        correct cell index based on these two numbers. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location from the cell with index \f$m\f$. In this class,
        the function first determines the radial and angular bin indices \f$i\f$ and \f$k\f$ that
        correspond to the cell index \f$m\f$. The spherical coordinates of the central position are
        subsequently determined from \f[ \begin{split} r &= \frac{r_i + r_{i+1}}{2} \\ \theta &=
        \frac{\theta_k + \theta_{k+1}}{2} \\ \phi &= 0. \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. In this class,
        the function first determines the radial and angular bin indices \f$i\f$ and \f$k\f$ that
        correspond to the cell index \f$m\f$. Then a random radius \f$r\f$, a random inclination
        \f$\theta\f$, and a random azimuth \f$\phi\f$ are determined using \f[ \begin{split} r &=
        \left( r_i^3 + {\cal{X}}_1\,(r_{i+1}^3-r_i^3) \right)^{1/3} \\ \cos\theta &= \cos\theta_k +
        {\cal{X}}_2\, (\cos\theta_{k+1}-\cos\theta_k) \\ \phi &= 2\pi\,{\cal{X}}_3, \end{split} \f]
        with \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a 2D spherical grid, implemented as a
        private PathSegmentGenerator subclass. The algorithm used to construct the path is
        described below.

        We represent the path by its parameter equation \f${\bf{x}}={\bf{r}}+s\,{\bf{k}}\f$, and we
        assume that \f${\bf{k}}\f$ is a unit vector. The two intersection points with a radial
        boundary sphere \f${\bf{x}}^2=R^2\f$ are obtained by solving the quadratic equation \f$s^2
        + 2\,({\bf{r}}\cdot{\bf{k}})\,s + ({\bf{r}}^2-R^2)=0\f$ for \f$s\f$. The two intersection
        points with an angular boundary cone \f$x_z^2=c^2\,{\bf{x}}^2\f$ (with \f$c=\cos\theta\f$)
        are obtained by solving the quadratic equation \f$(c^2-k_z^2)\,s^2 +
        2\,(c^2\,{\bf{r}}\cdot{\bf{k}}-r_z k_z)\,s + (c^2\,{\bf{r}}^2-r_z^2)=0\f$ for \f$s\f$. The
        intersection points with the reflected cone are always more distant than the other cell
        boundaries (the requirement to include the xy-plane \f$\theta=\pi/2\f$ in the grid ensures
        that this is true) and thus these phantom points are automatically ignored. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the xz plane to the specified
        SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

private:
    /** This private function returns the cell index corresponding to the radial index \f$i\f$ and
        the angular index \f$k\f$. The correspondence is simply \f$m=k+N_\theta\,i\f$. */
    int index(int i, int k) const;

    /** This private function calculates the radial index \f$i\f$ and the angular index \f$k\f$
        from a cell index \f$m\f$. As the correspondence between \f$m\f$, \f$i\f$ and \f$k\f$ is
        given by \f$m=k+N_\theta\,i\f$, one directly obtains \f$i=\lfloor m/N_\theta \rfloor\f$ and
        \f$k=m\!\mod N_\theta\f$. */
    void invertIndex(int m, int& i, int& k) const;

    //======================== Data Members ========================

private:
    int _Nr{0};
    int _Ntheta{0};
    Array _rv;
    Array _thetav;
    Array _cv;

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif

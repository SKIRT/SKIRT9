/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERE2DSPATIALGRID_HPP
#define SPHERE2DSPATIALGRID_HPP

#include "Mesh.hpp"
#include "SphereSpatialGrid.hpp"
class Random;

//////////////////////////////////////////////////////////////////////

/** The Sphere2DSpatialGrid class is subclass of the SphereSpatialGrid class, and represents
    axisymmetric spatial grids based on spherical coordinates. The grid is defined in the
    meridional plane and rotated around the z-axis. The meridional grid is specified through two
    sets of grid points:

    - \f$N_r+1\f$ radial grid points \f$r_i, \,i=0,\ldots,N_r\f$, with \f$0\le r_\text{min} =
    r_0\f$, \f$r_i<r_{i+1}\f$, and \f$r_{N_r} = r_\text{max}\f$.

    - \f$N_\theta+1\f$ angular (polar-inclination) grid points \f$\theta_j,
    \,j=0,\ldots,N_\theta\f$ with \f$0=\theta_0\f$, \f$\theta_j<\theta_{j+1}\f$, and
    \f$\theta_{N_\theta}=\pi\f$.

    \note The algorithm used by the path segment generator in this class requires that the xy-plane
    \f$\theta=\pi/2\f$ is included in the polar grid. If this is not the case, this point is
    automatically added, increasing the number of polar bins by one.

    In total there are \f$N_{\text{cells}} = N_r\,N_\theta\f$ cells in the grid. */
class Sphere2DSpatialGrid : public SphereSpatialGrid
{
    ITEM_CONCRETE(Sphere2DSpatialGrid, SphereSpatialGrid, "a 2D axisymmetric spatial grid in spherical coordinates")
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
        inclination angle cosines \f$c_j = \cos\theta_j\f$ for the boundary cones in the angular
        grid, to help speed up the path segment generator. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 2. */
    int dimension() const override;

    /** This function returns the number of cells \f$N_r\,N_\theta\f$ in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. It determines the radial
        and angular bin indices \f$i\f$ and \f$j\f$ corresponding to the cell index \f$m\f$, and
        then calculates the volume as \f[ V = \frac{2\pi}{3} \left(r_{i+1}^3-r_i^3\right)
        \left(\cos\theta_j-\cos\theta_{j+1}\right). \f] */
    double volume(int m) const override;

    /** This function returns the "diagonal" of the cell with index \f$m\f$. It determines the
        radial and angular bin indices \f$i\f$ and \f$j\f$ corresponding to the cell index \f$m\f$,
        and then calculates the distance between the outer/upper and inner/lower corners of the
        cell in the meridional plane. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. It determines the radial and angular bin indices \f$i\f$ and \f$j\f$ of the
        cell containing the position and calculates the correct cell index based on these two
        numbers. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. It determines
        the radial and angular bin indices \f$i\f$ and \f$j\f$ corresponding to the cell index
        \f$m\f$. The spherical coordinates of the central position are subsequently determined from
        \f[ \begin{split} r &= \frac{r_i + r_{i+1}}{2} \\ \theta &= \frac{\theta_j +
        \theta_{j+1}}{2} \\ \varphi &= 0. \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location in the cell with index \f$m\f$. It determines the
        radial and angular bin indices \f$i\f$ and \f$j\f$ corresponding to the cell index \f$m\f$.
        Then a random radius \f$r\f$, a random inclination \f$\theta\f$, and a random azimuth
        \f$\varphi\f$ are determined using \f[ \begin{split} r &= \left( r_i^3 +
        {\cal{X}}_1\,(r_{i+1}^3-r_i^3) \right)^{1/3} \\ \cos\theta &= \cos\theta_j + {\cal{X}}_2\,
        (\cos\theta_{j+1}-\cos\theta_j) \\ \varphi &= 2\pi\,{\cal{X}}_3, \end{split} \f] with
        \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a 2D spherical grid, implemented as a
        private PathSegmentGenerator subclass.

        We represent the path by its parameter equation \f${\bf{x}}={\bf{r}}+s\,{\bf{k}}\f$, and we
        assume that \f${\bf{k}}\f$ is a unit vector. The two intersection points with a radial
        boundary sphere \f${\bf{x}}^2=r_*^2\f$ are obtained by solving the quadratic equation
        \f$s^2 + 2\,({\bf{r}}\cdot{\bf{k}})\,s + ({\bf{r}}^2-r_*^2)=0\f$ for \f$s\f$. The two
        intersection points with an angular boundary cone \f$x_z^2=c^2\,{\bf{x}}^2\f$ (with
        \f$c=\cos\theta_*\f$) are obtained by solving the quadratic equation \f$(c^2-k_z^2)\,s^2 +
        2\,(c^2\,{\bf{r}}\cdot{\bf{k}}-r_z k_z)\,s + (c^2\,{\bf{r}}^2-r_z^2)=0\f$ for \f$s\f$. The
        intersection points with the reflected cone are always more distant than the other cell
        boundaries (the requirement to include the xy-plane \f$\theta=\pi/2\f$ in the grid ensures
        that this is true) and thus these phantom points are automatically ignored.

        The segment generator progresses the starting point of the path through the grid along the
        path's direction. For each step along the way, it calculates the distances to the
        intersections with all candidate borders of the current cell, and then selects the nearest
        intersection point. */
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
        the angular index \f$j\f$. The correspondence is simply \f$m=j+N_\theta\,i\f$. */
    int index(int i, int j) const;

    /** This function obtains the coordinates in the meridional plane for the corners of the cell
        with index \f$m\f$. It determines the radial and angular bin indices \f$i\f$ and \f$j\f$
        corresponding to the cell index \f$m\f$ using \f$i=\lfloor m/N_z \rfloor\f$ and
        \f$j=m\,\text{mod}\,N_z\f$.

        If both the resulting bin indices are within range, the function stores the corresponding
        cell corner coordinates in the provided arguments and returns true. If any of the indices
        are out of range, the function returns false and the contents of the provided arguments
        remains unchanged. */
    bool getCoords(int m, double& rmin, double& thetamin, double& rmax, double& thetamax) const;

    //======================== Data Members ========================

private:
    int _Nr{0};
    int _Ntheta{0};
    int _Ncells{0};
    Array _rv;
    Array _thetav;
    Array _cv;

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif

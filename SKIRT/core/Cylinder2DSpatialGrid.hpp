/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDER2DSPATIALGRID_HPP
#define CYLINDER2DSPATIALGRID_HPP

#include "Array.hpp"
#include "CylinderSpatialGrid.hpp"
#include "Mesh.hpp"

////////////////////////////////////////////////////////////////////

/** The Cylinder2DSpatialGrid class is subclass of the CylinderSpatialGrid class, and represents
    axisymmetric spatial grids based on cylindrical coordinates. The grid is defined in the
    meridional plane and rotated around the z-axis. The meridional grid is specified through two
    sets of grid points:

    - \f$N_R+1\f$ radial grid points \f$R_i, \,i=0,\ldots,N_R\f$, with \f$0\le R_\text{min} =
    R_0\f$, \f$R_i<R_{i+1}\f$, and \f$R_{N_R} = R_\text{max}\f$.

    - \f$N_z+1\f$ vertical grid points \f$z_k, \,k=0,\ldots,N_z\f$, with \f$z_\text{min} = z_0\f$,
    \f$z_i<z_{i+1}\f$, and \f$z_{N_z} = z_\text{max}\f$.

    In total there are \f$N_{\text{cells}} = N_R\,N_z\f$ cells in the grid. */
class Cylinder2DSpatialGrid : public CylinderSpatialGrid
{
    ITEM_CONCRETE(Cylinder2DSpatialGrid, CylinderSpatialGrid,
                  "a 2D axisymmetric spatial grid in cylindrical coordinates")
        ATTRIBUTE_TYPE_ALLOWED_IF(Cylinder2DSpatialGrid, "!Dimension3")

        PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

        PROPERTY_ITEM(meshZ, Mesh, "the bin distribution in the Z direction")
        ATTRIBUTE_DEFAULT_VALUE(meshZ, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up a number of data members that depend on the Mesh objects configured
        for this grid. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 2. */
    int dimension() const override;

    /** This function returns the number of cells \f$N_R\,N_z\f$ in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. It determines the radial
        and vertical bin indices \f$i\f$ and \f$k\f$ corresponding to the index \f$m\f$, and then
        calculates the volume as \f[ V = \pi (R_{i+1}^2-R_i^2) (z_{k+1}-z_k). \f] */
    double volume(int m) const override;

    /** This function returns the "diagonal" of the cell with index \f$m\f$. It determines the
        radial and vertical bin indices \f$i\f$ and \f$k\f$ corresponding to the index \f$m\f$, and
        then calculates the distance between the outer/upper and inner/lower corners of the cell in
        the meridional plane. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. It determines the radial and vertical bin indices \f$i\f$ and \f$k\f$ of
        the cell containing the position and calculates the correct cell index based on these two
        numbers. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. It determines
        the radial and vertical bin indices \f$i\f$ and \f$k\f$ corresponding to the index \f$m\f$.
        The cylindrical coordinates of the central position are subsequently determined from \f[
        \begin{split} R &= \frac{R_i + R_{i+1}}{2} \\ \phi &= 0 \\ z &= \frac{z_k + z_{k+1}}{2}
        \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location in the cell with index \f$m\f$. It determines the
        radial and vertical bin indices \f$i\f$ and \f$k\f$ corresponding to the index \f$m\f$, and
        then generates the cylindrical coordinates of a random position using \f[ \begin{split} R
        &= \sqrt{R_i^2 + {\cal{X}}_1\,(R_{i+1}^2-R_i^2)} \\ \phi &= 2\pi\,{\cal{X}}_2 \\ z &= z_k +
        {\cal{X}}_3\, (z_{k+1}-z_k), \end{split} \f] with \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and
        \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a 2D cylindrical grid, implemented as a
        private PathSegmentGenerator subclass.

        \image html Cylinder2DSpatialGrid.png

        The figure above is a projection on the equatorial plane. It shows coordinate system origin
        \f$\bf{o}\f$ and the path under consideration, defined by its projected starting point
        \f$\bf{R}\f$ and projected direction vector \f$\bf{k}_\text{q}\f$. Now consider the point
        of closest approach \f$\bf{c}\f$. The distance \f$q\f$ from \f$\bf{c}\f$ to \f$\bf{R}\f$
        can be obtained by projecting the position vector \f${\bf{R}}\f$ onto the line formed by
        the path. This distance is given by \f$q=\bf{R}\cdot\bf{k}_\text{q}/||\bf{k}_\text{q}||\f$,
        where we divided the dot product by the norm of the projected direction vector because it
        is not normalized. The resulting value is negative if \f$\bf{R}\f$ is before \f$\bf{c}\f$
        (i.e. it is going inward) and positive if \f$\bf{R}\f$ is after \f$\bf{c}\f$ (i.e. it is
        going outward). We further define the impact parameter \f$p\f$ as the distance of closest
        approach. From the rectangular triangles illustrated in the figure it is easily seen that
        \f$p^2 + q^2 = R^2\f$ and \f$p^2 + q_*^2 = R_*^2\f$. To obtain the actual 3D distance
        traveled, the \f$q\f$ values must be divided by the projected direction cosine
        \f$||\bf{k}_\text{q}||\f$.

        In the vertical direction, the distance to the intersection of the path with a horizontal
        plane at height \f$z_*\f$ is given by \f$(z_*-z)/k_\text{z}\f$ as usual.

        The segment generator uses these relations to step the path through the cells, finding the
        nearest intersection and updating the values of \f$q\f$, \f$\bf{R}\f$ and \f$z\f$ along the
        way as a proxy for updating the 3D position \f${\bf{r}}\f$. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the xz plane to the specified
        SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

private:
    /** This function returns the index \f$m\f$ corresponding to the radial index \f$i\f$ and the
        vertical index \f$k\f$. The correspondence is simply \f$m=k+N_z\,i\f$. */
    int index(int i, int k) const;

    /** This function obtains the coordinates in the meridional plane for the corners of the cell
        with index \f$m\f$. It determines the radial and vertical bin indices \f$i\f$ and \f$k\f$
        corresponding to the index \f$m\f$ using \f$i=\lfloor m/N_z \rfloor\f$ and
        \f$k=m\,\text{mod}\,N_z\f$.

        If both the resulting bin indices are within range, the function stores the corresponding
        cell corner coordinates in the provided arguments and returns true. If any of the indices
        are out of range, the function returns false and the contents of the provided arguments
        remains unchanged. */
    bool getCoords(int m, double& Rmin, double& zmin, double& Rmax, double& zmax) const;

    //======================== Data Members ========================

private:
    int _NR{0};
    int _Nz{0};
    int _Ncells{0};
    Array _Rv;
    Array _zv;

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif

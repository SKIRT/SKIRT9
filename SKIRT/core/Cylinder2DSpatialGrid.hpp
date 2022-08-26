/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDER2DSPATIALGRID_HPP
#define CYLINDER2DSPATIALGRID_HPP

#include "Array.hpp"
#include "CylinderSpatialGrid.hpp"
#include "MoveableMesh.hpp"

////////////////////////////////////////////////////////////////////

/** The Cylinder2DSpatialGrid class is subclass of the CylinderSpatialGrid class, and represents
    axisymmetric spatial grids based on cylindrical coordinates. The grid is defined in the
    meridional plane and rotated around the Z-axis. The meridional grid is specified through a set
    of \f$N_R+1\f$ radial grid points \f$R_i\f$ (with \f$i=0,\ldots,N_R\f$) and a set of
    \f$N_z+1\f$ vertical grid points \f$z_k\f$ (with \f$k=0,\ldots,N_z\f$). In total there are
    \f$N_{\text{cells}} = N_R\,N_z\f$ cells in the Spatial grid. */
class Cylinder2DSpatialGrid : public CylinderSpatialGrid
{
    ITEM_CONCRETE(Cylinder2DSpatialGrid, CylinderSpatialGrid, "an axisymmetric spatial grid in cylindrical coordinates")
        ATTRIBUTE_TYPE_ALLOWED_IF(Cylinder2DSpatialGrid, "!Dimension3")

        PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

        PROPERTY_ITEM(meshZ, MoveableMesh, "the bin distribution in the Z direction")
        ATTRIBUTE_DEFAULT_VALUE(meshZ, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up a number of data members that depend on the Mesh objects configured
        for this grid. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 2 for this class. */
    int dimension() const override;

    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. For an axisymmetric grid,
        the function determines the radial and vertical bin indices \f$i\f$ and \f$k\f$ that
        correspond to the index \f$m\f$, and then calculates the volume as \f[ V = \pi
        \left(R_{i+1}-R_i\right)^2 \left(z_{k+1}-z_k\right). \f] */
    double volume(int m) const override;

    /** This function returns the "diagonal" of the cell with index \f$m\f$. For 2D cylindrical
        grids, it returns the distance between the outer/upper and inner/lower corners of the cell
        in the meridional plane. */
    double diagonal(int m) const override;

    /** This function returns the number of the cell that contains the position \f${\bf{r}}\f$. It
        just determines the radial and vertical bin indices and calculates the correct cell index
        based on these two numbers. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location from the cell with index \f$m\f$. For an
        axisymmetric grid, the function first determines the radial and vertical bin indices
        \f$i\f$ and \f$k\f$ that correspond to the index \f$m\f$. The cylindrical coordinates of
        the central position are subsequently determined from \f[ \begin{split} R &= \frac{R_i +
        R_{i+1}}{2} \\ \phi &= 0 \\ z &= \frac{z_k + z_{k+1}}{2} \end{split} \f] A position with
        these cylindrical coordinates is returned. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. For an
        axisymmetric grid, the function first determines the radial and vertical bin indices
        \f$i\f$ and \f$k\f$ that correspond to the index \f$m\f$. Then a random radius \f$R\f$, a
        random azimuth \f$\phi\f$, and a random height \f$z\f$ are determined using \f[
        \begin{split} R &= \sqrt{R_i^2 + {\cal{X}}_1\,(R_{i+1}-R_i)^2} \\ \phi &= 2\pi\,{\cal{X}}_2
        \\ z &= z_k + {\cal{X}}_3\, (z_{k+1}-z_k), \end{split} \f] with \f${\cal{X}}_1\f$,
        \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three uniform deviates. A position with these
        cylindrical coordinates is returned. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for an axisymmetric grid, implemented as a
        private PathSegmentGenerator subclass.

        The algorithm used to construct the path is fairly straightforward because of the symmetry
        in the grid. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the xz plane to the specified
        SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

private:
    /** This private function returns the index corresponding to the radial index \f$i\f$ and the
        vertical index \f$k\f$. The correspondence is simply \f$m=k+N_z\,i\f$. */
    int index(int i, int k) const;

    /** This private function calculates the radial index \f$i\f$ and the vertical index \f$k\f$
        from a cell index \f$m\f$. As the correspondence between \f$m\f$, \f$i\f$ and \f$k\f$ is
        given by \f$m=k+N_z\,i\f$, one directly obtains \f$i=\lfloor m/N_z \rfloor\f$ and
        \f$k=m\!\mod N_z\f$. */
    void invertIndex(int m, int& i, int& k) const;

    //======================== Data Members ========================

private:
    int _NR{0};
    int _Nz{0};
    Array _Rv;
    Array _zv;

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif

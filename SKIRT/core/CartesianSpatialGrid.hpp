/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CARTESIANSPATIALGRID_HPP
#define CARTESIANSPATIALGRID_HPP

#include "Array.hpp"
#include "BoxSpatialGrid.hpp"
#include "MoveableMesh.hpp"

////////////////////////////////////////////////////////////////////

/** The CartesianSpatialGrid class is subclass of the BoxSpatialGrid class, and represents
    three-dimensional spatial grids based on a regular Cartesian grid. Each cell in such a grid is
    a little cuboid (not necessarily all with the same size or axis ratios). */
class CartesianSpatialGrid : public BoxSpatialGrid
{
    ITEM_CONCRETE(CartesianSpatialGrid, BoxSpatialGrid, "a Cartesian spatial grid")

        PROPERTY_ITEM(meshX, MoveableMesh, "the bin distribution in the X direction")
        ATTRIBUTE_DEFAULT_VALUE(meshX, "LinMesh")

        PROPERTY_ITEM(meshY, MoveableMesh, "the bin distribution in the Y direction")
        ATTRIBUTE_DEFAULT_VALUE(meshY, "LinMesh")

        PROPERTY_ITEM(meshZ, MoveableMesh, "the bin distribution in the Z direction")
        ATTRIBUTE_DEFAULT_VALUE(meshZ, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets up a number of data members that depend on the Mesh objects configured
        for this grid. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. For a cartesian grid, the
        function determines the bin indices \f$i\f$, \f$j\f$ and \f$k\f$ corresponding to the X, Y
        and Z directions. The volume is then easily calculated as \f$V = (x_{i+1}-x_i)\,
        (y_{j+1}-y_j)\, (z_{k+1}-z_k) \f$. */
    double volume(int m) const override;

    /** This function returns the diagonal of the cell with index \f$m\f$. For a cartesian grid, the
        function determines the bin indices \f$i\f$, \f$j\f$ and \f$k\f$ corresponding to the X, Y
        and Z directions. The diagonal is then easily calculated as \f$d = \sqrt{ (x_{i+1}-x_i)^2 +
        (y_{j+1}-y_j)^2 + (z_{k+1}-z_k)^2 }\f$. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. For a cartesian grid, the function determines the bin indices in the X, Y
        and Z directions and calculates the correct index based on these indices. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location from the cell with index \f$m\f$. For a
        cartesian grid, the function first determines the bin indices \f$i\f$ and \f$k\f$ in the X,
        Y and Z directions that correspond to the index \f$m\f$. Then the coordinates \f$x\f$,
        \f$y\f$ and \f$z\f$ are determined using \f[ \begin{split} x &= \frac{x_i + x_{i+1}}{2} \\
        y &= \frac{y_j + y_{j+1}}{2} \\ z &= \frac{z_k + z_{k+1}}{2} \end{split} \f] A position
        with these cartesian coordinates is returned. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. For a cartesian
        grid, the function first determines the bin indices \f$i\f$ and \f$k\f$ in the X, Y and Z
        directions that correspond to the index \f$m\f$. Then random coordinates \f$x\f$, \f$y\f$
        and \f$z\f$ are determined using \f[ \begin{split} x &= x_i + {\cal{X}}_1\,(x_{i+1}-x_i) \\
        y &= y_j + {\cal{X}}_2\,(y_{j+1}-y_j) \\ z &= z_k + {\cal{X}}_3\, (z_{k+1}-z_k),
        \end{split} \f] with \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three
        uniform deviates. A position with these cartesian coordinates is returned. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a cartesian grid, implemented as a private
        PathSegmentGenerator subclass.

        The algorithm used to construct the path is fairly straightforward because all cells are
        cuboids lined up with the coordinate axes and the neighboring cells are easily found by
        manipulating cell indices. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid structure with the xy plane to the
        specified SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid structure with the xz plane to the
        specified SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid structure with the yz plane to the
        specified SpatialGridPlotFile object. */
    void write_yz(SpatialGridPlotFile* outfile) const override;

    /** This function writes 3D information for all cells in the grid structure to the specified
        SpatialGridPlotFile object. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

private:
    /** This function returns the index \f$m\f$ corresponding to the three bin indices \f$i\f$,
       \f$j\f$ and \f$k\f$. The correspondence is \f$m=k+j\,N_z+i\,N_y\,N_z\f$. */
    int index(int i, int j, int k) const;

    /** This function calculates the three bin indices \f$i\f$, \f$j\f$ and \f$k\f$ of the cell
        index \f$m\f$, and then returns the coordinates of the corresponding cell as a Box object.
        Since the relation between the index and the three bin indices is
        \f$m=k+j\,N_z+i\,N_y\,N_z\f$, we can use the formulae \f[ \begin{split} i &= \lfloor
        m/(N_y\,N_z) \rfloor \\ j &= \lfloor (m-i\,N_y\,N_z)/N_z \rfloor \\ k &=
        m\,{\text{mod}}\,N_z. \end{split} \f] */
    Box box(int m) const;

    //======================== Data Members ========================

private:
    int _Nx{0};
    int _Ny{0};
    int _Nz{0};
    Array _xv;
    Array _yv;
    Array _zv;

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

////////////////////////////////////////////////////////////////////

#endif

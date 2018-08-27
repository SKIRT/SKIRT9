/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHSPATIALGRID_HPP
#define ADAPTIVEMESHSPATIALGRID_HPP

#include "SpatialGrid.hpp"
#include "DensityInCellInterface.hpp"
class AdaptiveMeshSnapshot;

////////////////////////////////////////////////////////////////////

/** An instance of the AdaptiveMeshSpatialGrid class represents a three-dimensional spatial grid,
    the structure of which is described by an imported adaptive mesh. In fact, this class directly
    uses the adaptive mesh created by an AdaptiveMeshGeometry or AdaptiveMeshMedium object. For
    this to work, the medium system must include at least one component that uses an imported
    AdaptiveMeshSnapshot to represent its spatial distribution. */
class AdaptiveMeshSpatialGrid : public SpatialGrid, public DensityInCellInterface
{
    ITEM_CONCRETE(AdaptiveMeshSpatialGrid, SpatialGrid, "a spatial grid based on an imported adaptive mesh snapshot")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function locates the adaptive mesh in the simulation hierarchy, and remembers a
        pointer to it. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is always 3 for an adaptive mesh. */
    int dimension() const override;

    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the bounding box that encloses the grid. */
    Box boundingBox() const override;

    /** This function returns the volume of the cell with index \f$m\f$. */
    double volume(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. */
    Position randomPositionInCell(int m) const override;

    /** This function implements the DensityInCellInterface interface. It returns the density for
        medium component \f$h\f$ in the grid cell with index \f$m\f$. For an adaptive mesh grid,
        this interface can be offered only when the medium system consists of a single component,
        which then, by definition, supplies the adaptive mesh. Thus, the component index \f$h\f$
        passed to this function should always be zero; in fact, its value is actually ignored. */
    double density(int h, int m) const override;

    // TO DO: allow DensityInCellInterface only when there is just a single medium component

    /** This function calculates a path through the grid. The SpatialGridPath object passed as an
        argument specifies the starting position \f${\bf{r}}\f$ and the direction \f${\bf{k}}\f$
        for the path. The data on the calculated path are added back into the same object. */
    void path(SpatialGridPath* path) const override;

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

    //======================== Data Members ========================

private:
    AdaptiveMeshSnapshot* _snapshot{nullptr};  // adaptive grid obtained from medium component
    double _norm{0.};                          // normalization factor obtained from medium component
};

////////////////////////////////////////////////////////////////////

#endif

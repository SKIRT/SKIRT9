/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHSPATIALGRID_HPP
#define ADAPTIVEMESHSPATIALGRID_HPP

#include "DensityInCellInterface.hpp"
#include "SpatialGrid.hpp"
class AdaptiveMeshSnapshot;

////////////////////////////////////////////////////////////////////

/** An instance of the AdaptiveMeshSpatialGrid class represents a three-dimensional spatial grid,
    the structure of which is described by an imported adaptive mesh. In fact, this class directly
    uses the adaptive mesh created by an AdaptiveMeshGeometry or AdaptiveMeshMedium object. For
    this to work, the medium system must include at least one component that uses an imported
    AdaptiveMeshSnapshot to represent its spatial distribution. If multiple media components are
    based on an imported AdaptiveMeshSnapshot, the first component (in configuration order) is used
    for the spatial grid. */
class AdaptiveMeshSpatialGrid : public SpatialGrid, public DensityInCellInterface
{
    ITEM_CONCRETE(AdaptiveMeshSpatialGrid, SpatialGrid, "a spatial grid taken from an imported adaptive mesh snapshot")
        ATTRIBUTE_TYPE_ALLOWED_IF(AdaptiveMeshSpatialGrid, "AdaptiveMeshInterface")
        ATTRIBUTE_TYPE_DISPLAYED_IF(AdaptiveMeshSpatialGrid, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function locates the adaptive mesh snapshot in the simulation hierarchy, and remembers
        a pointer to it. */
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

    /** This function returns the diagonal of the cell with index \f$m\f$. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this spatial grid type. For the adaptive
        mesh grid, the path segment generator is actually implemented in the AdaptiveMeshSnapshot
        class. */
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

public:
    /** This function implements the DensityInCellInterface interface. It returns the number
        density for medium component \f$h\f$ in the grid cell with index \f$m\f$. For an adaptive
        mesh grid, this interface can be offered only when the medium system consists of a single
        component, which then, by definition, supplies the adaptive mesh. Thus, the component index
        \f$h\f$ passed to this function should always be zero; in fact, its value is actually
        ignored. */
    double numberDensity(int h, int m) const override;

protected:
    /** This function is used by the interface() function to ensure that the receiving item can
        actually offer the specified interface. If the requested interface is the
        DensityInCellInterface, the implementation in this class returns true if the medium system
        consists of a single component, and false otherwise. For other requested interfaces, the
        function invokes its counterpart in the base class. */
    bool offersInterface(const std::type_info& interfaceTypeInfo) const override;

    //======================== Data Members ========================

private:
    AdaptiveMeshSnapshot* _mesh{nullptr};  // adaptive mesh snapshot obtained from medium component
    double _norm{0.};                      // normalization factor obtained from medium component
};

////////////////////////////////////////////////////////////////////

#endif

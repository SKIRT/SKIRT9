/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIMESHSPATIALGRID_HPP
#define VORONOIMESHSPATIALGRID_HPP

#include "BoxSpatialGrid.hpp"
#include "DensityInCellInterface.hpp"
class VoronoiMeshSnapshot;

//////////////////////////////////////////////////////////////////////

/** VoronoiMeshSpatialGrid is a concrete subclass of the SpatialGrid class. It represents a
    three-dimensional grid based on a Voronoi tesselation of the cuboidal spatial domain of the
    simulation. See the VoronoiMeshSnapshot class for more information on Voronoi tesselations.

    The class offers several options for determining the positions of the sites generating the
    Voronoi tesselation. A specified number of sites can be distributed randomly over the domain,
    either uniformly or with the same overall density distribution as the medium. Alternatively,
    the positions can be copied from the sites in the imported distribution(s).

    Furthermore, the user can opt to perform a relaxation step on the site positions to avoid
    overly elongated cells. */
class VoronoiMeshSpatialGrid : public BoxSpatialGrid, public DensityInCellInterface
{
    /** The enumeration type indicating the policy for determining the positions of the sites. */
    ENUM_DEF(Policy, Uniform, CentralPeak, DustDensity, ElectronDensity, GasDensity, File, ImportedSites, ImportedMesh)
        ENUM_VAL(Policy, Uniform, "random from uniform distribution")
        ENUM_VAL(Policy, CentralPeak, "random from distribution with a steep central peak")
        ENUM_VAL(Policy, DustDensity, "random from dust density distribution")
        ENUM_VAL(Policy, ElectronDensity, "random from electron density distribution")
        ENUM_VAL(Policy, GasDensity, "random from gas density distribution")
        ENUM_VAL(Policy, File, "loaded from text column data file")
        ENUM_VAL(Policy, ImportedSites, "positions of particles, sites or cells in imported distribution")
        ENUM_VAL(Policy, ImportedMesh, "employ imported Voronoi mesh in medium system")
    ENUM_END()

    ITEM_CONCRETE(VoronoiMeshSpatialGrid, BoxSpatialGrid, "a Voronoi tessellation-based spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(VoronoiMeshSpatialGrid, "Level2")

        PROPERTY_ENUM(policy, Policy, "the policy for determining the positions of the sites")
        ATTRIBUTE_DEFAULT_VALUE(policy, "DustDensity")

        PROPERTY_INT(numSites, "the number of random sites (or cells in the grid)")
        ATTRIBUTE_MIN_VALUE(numSites, "5")
        ATTRIBUTE_DEFAULT_VALUE(numSites, "500")
        ATTRIBUTE_RELEVANT_IF(numSites, "policyUniform|policyCentralPeak|policyDustDensity|"
                                        "policyElectronDensity|policyGasDensity")

        PROPERTY_STRING(filename, "the name of the file containing the site positions")
        ATTRIBUTE_RELEVANT_IF(filename, "policyFile")

        PROPERTY_BOOL(relaxSites, "perform site relaxation to avoid overly elongated cells")
        ATTRIBUTE_DEFAULT_VALUE(relaxSites, "false")
        ATTRIBUTE_RELEVANT_IF(relaxSites, "!policyImportedMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor releases the Voronoi mesh if this object owns it. */
    ~VoronoiMeshSpatialGrid();

protected:
    /** This function verifies that the attributes have been appropriately set, generates or
        retrieves the site positions for constructing the Voronoi tessellation according to the
        configured policy, and finally constructs the Voronoi tessellation through an instance of
        the VoronoiMeshSnapshot class. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. */
    double volume(int m) const override;

    /** This function returns the approximate diagonal of the cell with index \f$m\f$. For a
        Voronoi grid, it returns \f$(3V)^(1/3)\f$ where \f$V\f$ is the volume of the cell. This
        corresponds to the correct diagonal only for cubical cells. */
    double diagonal(int m) const override;

    /** This function returns the index of the cell that contains the position \f${\bf{r}}\f$. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. In this class
        the function returns the centroid of the Voronoi cell. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this spatial grid type. For the Voronoi
        mesh grid, the path segment generator is actually implemented in the VoronoiMeshSnapshot
        class. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

    /** This function outputs the grid plot files; it is provided here because the regular
        mechanism does not apply. The function reconstructs the Voronoi tesselation in order to
        produce the coordinates of the Voronoi cell vertices. */
    void writeGridPlotFiles(const SimulationItem* probe) const override;

    /** This function implements the DensityInCellInterface interface. It returns the number
        density for medium component \f$h\f$ in the grid cell with index \f$m\f$. For a Voronoi
        mesh grid, this interface can be offered only if the "ImportedMesh" policy has been
        configured, and the medium system consists of a single component, which then, by
        definition, supplies the Voronoi mesh. Thus, the component index \f$h\f$ passed to this
        function should always be zero; in fact, its value is actually ignored. */
    double numberDensity(int h, int m) const override;

protected:
    /** This function is used by the interface() function to ensure that the receiving item can
        actually offer the specified interface. If the requested interface is the
        DensityInCellInterface, the implementation in this class returns true if the "ImportedMesh"
        policy has been configured and the medium system consists of a single component, and false
        otherwise. For other requested interfaces, the function invokes its counterpart in the base
        class. */
    bool offersInterface(const std::type_info& interfaceTypeInfo) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    VoronoiMeshSnapshot* _mesh{nullptr};  // Voronoi mesh snapshot created here or obtained from medium component
    double _norm{0.};                     // in the latter case, normalization factor obtained from medium component
};

//////////////////////////////////////////////////////////////////////

#endif

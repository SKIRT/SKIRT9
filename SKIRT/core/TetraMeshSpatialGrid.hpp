/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TETRAMESHSPATIALGRID_HPP
#define TETRAMESHSPATIALGRID_HPP

#include "BoxSpatialGrid.hpp"
#include "DensityInCellInterface.hpp"
#include "Medium.hpp"
class TetraMeshSnapshot;

//////////////////////////////////////////////////////////////////////

/** TetraMeshSpatialGrid is a concrete subclass of the SpatialGrid class. It represents a
    three-dimensional grid based on a Tetra tesselation of the cuboidal spatial domain of the
    simulation. See the TetraMeshSnapshot class for more information on Tetra tesselations.

    The class offers several options for determining the positions of the sites generating the
    Tetra tesselation. A specified number of sites can be distributed randomly over the domain,
    either uniformly or with the same overall density distribution as the medium. Alternatively,
    the positions can be copied from the sites in the imported distribution(s).

    Furthermore, the user can opt to perform a relaxation step on the site positions to avoid
    overly elongated cells. */
class TetraMeshSpatialGrid : public BoxSpatialGrid, public DensityInCellInterface
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
        ENUM_VAL(Policy, ImportedMesh, "employ imported Tetra mesh in medium system")
    ENUM_END()

    ITEM_CONCRETE(TetraMeshSpatialGrid, BoxSpatialGrid, "a Tetra tessellation-based spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TetraMeshSpatialGrid, "Level2")

        PROPERTY_ENUM(policy, Policy, "the policy for determining the positions of the sites")
        ATTRIBUTE_DEFAULT_VALUE(policy, "DustDensity")

        PROPERTY_DOUBLE(maxDustFraction, "the maximum fraction of dust contained in each cell")
        ATTRIBUTE_MIN_VALUE(maxDustFraction, "[0")
        ATTRIBUTE_MAX_VALUE(maxDustFraction, "1e-2]")
        ATTRIBUTE_DEFAULT_VALUE(maxDustFraction, "1e-6")
        ATTRIBUTE_DISPLAYED_IF(maxDustFraction, "DustMix")

        PROPERTY_DOUBLE(maxDustOpticalDepth, "the maximum diagonal dust optical depth for each cell")
        ATTRIBUTE_MIN_VALUE(maxDustOpticalDepth, "[0")
        ATTRIBUTE_MAX_VALUE(maxDustOpticalDepth, "100]")
        ATTRIBUTE_DEFAULT_VALUE(maxDustOpticalDepth, "0")
        ATTRIBUTE_DISPLAYED_IF(maxDustOpticalDepth, "DustMix&Level2")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to evaluate the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_RELEVANT_IF(wavelength, "maxDustOpticalDepth")

        PROPERTY_DOUBLE(maxDustDensityDispersion, "the maximum dust density dispersion in each cell")
        ATTRIBUTE_MIN_VALUE(maxDustDensityDispersion, "[0")
        ATTRIBUTE_MAX_VALUE(maxDustDensityDispersion, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxDustDensityDispersion, "0")
        ATTRIBUTE_DISPLAYED_IF(maxDustDensityDispersion, "DustMix&Level2")

        PROPERTY_DOUBLE(maxElectronFraction, "the maximum fraction of electrons contained in each cell")
        ATTRIBUTE_MIN_VALUE(maxElectronFraction, "[0")
        ATTRIBUTE_MAX_VALUE(maxElectronFraction, "1e-2]")
        ATTRIBUTE_DEFAULT_VALUE(maxElectronFraction, "1e-6")
        ATTRIBUTE_DISPLAYED_IF(maxElectronFraction, "ElectronMix")

        PROPERTY_DOUBLE(maxGasFraction, "the maximum fraction of gas contained in each cell")
        ATTRIBUTE_MIN_VALUE(maxGasFraction, "[0")
        ATTRIBUTE_MAX_VALUE(maxGasFraction, "1e-2]")
        ATTRIBUTE_DEFAULT_VALUE(maxGasFraction, "1e-6")
        ATTRIBUTE_DISPLAYED_IF(maxGasFraction, "GasMix")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor releases the Tetra mesh if this object owns it. */
    ~TetraMeshSpatialGrid();

protected:
    /** This function verifies that the attributes have been appropriately set, generates or
        retrieves the site positions for constructing the Tetra tessellation according to the
        configured policy, and finally constructs the Tetra tessellation through an instance of
        the TetraMeshSnapshot class. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    bool tetUnsuitable(double* pa, double* pb, double* pc, double* pd, double vol) const;

    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. */
    double volume(int m) const override;

    /** This function returns the approximate diagonal of the cell with index \f$m\f$. For a
        Tetra grid, it returns \f$(3V)^(1/3)\f$ where \f$V\f$ is the volume of the cell. This
        corresponds to the correct diagonal only for cubical cells. */
    double diagonal(int m) const override;

    /** This function returns the index of the cell that contains the position \f${\bf{r}}\f$. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. In this class
        the function returns the centroid of the Tetra cell. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this spatial grid type. For the Tetra
        mesh grid, the path segment generator is actually implemented in the TetraMeshSnapshot
        class. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

    /** This function outputs the grid plot files; it is provided here because the regular
        mechanism does not apply. The function reconstructs the Tetra tesselation in order to
        produce the coordinates of the Tetra cell vertices. */
    void writeGridPlotFiles(const SimulationItem* probe) const override;

    /** This function implements the DensityInCellInterface interface. It returns the number
        density for medium component \f$h\f$ in the grid cell with index \f$m\f$. For a Tetra
        mesh grid, this interface can be offered only if the "ImportedMesh" policy has been
        configured, and the medium system consists of a single component, which then, by
        definition, supplies the Tetra mesh. Thus, the component index \f$h\f$ passed to this
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
    TetraMeshSnapshot* _mesh{nullptr};  // Tetra mesh snapshot created here or obtained from medium component
    double _norm{0.};                   // in the latter case, normalization factor obtained from medium component

    // data members initialized by setupSelfBefore()
    Random* _random{nullptr};
    int _numSamples{0};

    // lists of medium components of each material type;
    // list remains empty if no criteria are enabled for the corresponding material type
    vector<Medium*> _dustMedia;
    vector<Medium*> _electronMedia;
    vector<Medium*> _gasMedia;

    // flags become true if corresponding criterion is enabled
    // (i.e. configured maximum is nonzero and material type is present)
    bool _hasAny{false};
    bool _hasDustAny{false};
    bool _hasDustFraction{false};
    bool _hasDustOpticalDepth{false};
    bool _hasDustDensityDispersion{false};
    bool _hasElectronFraction{false};
    bool _hasGasFraction{false};

    // cached values for each material type (valid if corresponding flag is enabled)
    double _dustMass{0.};
    double _dustKappa{0.};
    double _electronNumber{0.};
    double _gasNumber{0.};
};

//////////////////////////////////////////////////////////////////////

#endif

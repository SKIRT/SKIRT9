/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MOLECULARLINEGASMIX_HPP
#define MOLECULARLINEGASMIX_HPP

#include "EmittingGasMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MolecularLineGasMix class describes the material properties related to selected rotational,
    vibrational and/or rovibrational transitions in selected molecules and atoms. The class
    properties allow configuring the species and number of transitions to be considered.

    For each supported transition, the emission luminosity and absorption opacity in a given cell
    are determined from the gas properties defined in the input model and the local radiation field
    calculated by the simulation. The class performs an iterative, non-LTE calculation including
    the effects of collisional transitions (excitation and de-excitation with one or more types of
    interaction partners) and photonic transitions (spontaneous emission, absorption, and induced
    emission). This allows establishing the energy level distribution (population) for a wide range
    of material densities and radiation fields.

    <b>Supported species and transitions</b>

    The current implementation supports the following molecular or atomic species:

    - \c TwoLevelBenchmarkMolecule: a fictive test molecule with two energy levels and a single
    collisional partner as defined by van Zadelhoff et al. 2002 for the first benchmark problem
    described there. The corresponding line is at 1666.67 \f$\mu\mathrm{m}\f$.

    [TO DO]

    <b>Configuring the simulation</b>

    Simulations that include gas represented by the MolecularLineGasMix often also include dust,
    although this is not a requirement. In any case, the simulation must have one or more primary
    sources that trigger the molecular lines directly (e.g. the cosmic microwave background) or
    indirectly (e.g. by heating the dust and thus causing thermal dust emission), or both. The \em
    simulationMode should be set to "GasEmission" or "DustAndGasEmission" correspondingly.

    Because the secondary emission and the radiation field are calculated self-consistently, \em
    iterateSecondaryEmission should always be enabled. If the wavelength range of the primary
    source(s) significantly overlaps the emission lines of the species under consideration (in
    other words, if the opacity at the line wavelengths significantly affects the primary radiation
    field), then \em includePrimaryEmission in IterationOptions must be enabled as well. If not,
    this property can be left disabled. In both cases, \em iteratePrimaryEmission can be left
    disabled unless iteration is required for other media (for example, to self-consistently
    calculate radiative dust destruction).

    The radiation field wavelength grid should properly resolve the UV, optical, and infrared
    wavelength range (in case the simulation includes dust) and the wavelength ranges of the
    supported transition lines. Separate instruments can be configured for the relevant wavelength
    ranges, e.g. using a logarithmic grid for the continuum spectrum and linear grids for the line
    profiles.

    During the initial primary emission segment, the simulation determines the radiation field
    resulting from the primary sources and dust attenuation (molecular line absorption is taken to
    be zero at this stage). The resulting radiation field allows a first estimation of the level
    populations for each spatial grid cell; this calculation happens in the updateSpecificState()
    function.

    During secondary emission, the simulation takes into account emission from all configured media
    in addition to absorption by these same media, including molecular lines in both cases. The
    previously stored level populations allow calculating the line emission spectrum and the line
    absorption cross sections. This results in an updated radiation field, which will in turn
    influence the level populations (and for high optical depths, possibly the dust temperature),
    which in turn influences the secondary emission spectra. In order to obtain a self-consistent
    result, the simulation must therefore iterate over secondary emission.

    The input model must provide values for the spatial distribution of several medium properties,
    including the number density of the species under consideration, the number density of any
    relevant collisional partner species, the kinetic gas temperature, and the microturbulence
    level. These values remain constant during the simulation. Most often, this information will be
    read from an input file by associating the MolecularLineGasMix with a subclass of
    ImportedMedium. For that medium component, the ski file attribute \em importTemperature
    <b>must</b> be set to 'true', and \em importMetallicity and \em importVariableMixParams must be
    left at 'false'. The additional columns required by the material mix are automatically imported
    and are expected <b>after</b> all other columns. For example, if bulk velocities are also
    imported for this medium component (i.e. \em importVelocity is 'true'), the column order would
    be \f[ ..., n_\mathrm{mol}, T_\mathrm{kin}, v_\mathrm{x}, v_\mathrm{y}, v_\mathrm{z},
    n_\mathrm{col1} [, n_\mathrm{col1}, ...], v_\mathrm{turb}\f]

    For basic testing purposes, the MolecularLineGasMix can also be associated with a geometric
    medium component. The geometry then defines the spatial density distribution of the species
    under consideration (i.e. \f$n_\mathrm{mol}\f$), and the MolecularLineGasMix configuration
    properties specify a fixed default value for the other properties that will be used across the
    spatial domain. In this case, the number densities of the collisional partners are defined by a
    constant multiplier relative to \f$n_\mathrm{mol}\f$.

    <b>Level populations</b>

    The level populations \f$n_i\f$ (used later to calculate the emission luminosities and the
    absorption opacities) are obtained by solving the statistical equilibrium equations. The
    equation for the energy level with index \f$i\f$ (where indices increase from lowest to highest
    energy state) is given by

    \f[ \sum^{N_\mathrm{lp}}_{j>i} n_jA_{ji} + \sum^{N_\mathrm{lp}}_{j\neq i} \Big[ (n_jB_{ji} -
    n_iB_{ij})J_\lambda(\lambda_{ij})\Big] + \sum^{N_\mathrm{lp}}_{j\neq i}
    \big[n_jC_{ji}(n_\mathrm{col} ,T_\mathrm{kin}) - n_iC_{ij}(n_\mathrm{col}
    ,T_\mathrm{kin})\big]=0, \f]

    where \f$\lambda_{ij}\f$ is the transition wavelength and \f$A_{ij}\f$, \f$B_{ij}\f$,
    \f$B_{ji}\f$ are the Einstein coefficients for spontaneous emission, induced emission, and
    absorption, and \f$C_{ij}\f$, \f$C_{ji}\f$ the coefficients for collisional excitation and
    de-excitation. These coefficients are taken from the literature (references to be provided).

    \em Note

    We assume that the level populations of all molecules in a given spatial cell are the same, and
    to calculate them we use the mean intensity of the radiation field at the rest wavelength of
    each transition. This means that the radiation field wavelength grid needs just a single bin
    for each transition. We consider this assumption to be valid if the width of the velocity bin
    (or wavelength bin) of the observed spectrum is larger than the velocity dispersion of each
    cell. In reality, however, the molecules in a cell have a velocity dispersion and the perceived
    transition wavelength is different for each molecule. Therefore, one should calculate the mean
    intensities of each level transition at several velocity bins in the wavelength range
    corresponding to the velocity dispersion. However, the computational cost becomes very large.

    <b>Emission</b>

    The integrated line luminosity \f$L_\ell\f$ corresponding to the transition with index
    \f$\ell\f$ for a given spatial cell is given by

    \f[ L_\ell = \frac{hc}{\lambda_\ell} A_{ul} n_u V_\mathrm{cell}, \f]

    where \f$\lambda_\ell\f$ is the transition wavelength, \f$u\f$ and \f$l\f$ are the indices of
    the energy levels before and after the corresponding transition, and \f$V_\mathrm{cell}\f$ is
    the volume of the cell. The SKIRT framework automatically adds a Gaussian line profile assuming
    a thermal velocity corresponding to the kinetic temperature in the cell and the mass of a
    molecule, in addition to any Doppler shifts caused by the bulk velocity in the cell.

    <b>Absorption</b>

    During primary emission, absorption cannot be calculated (and is taken to be zero) because the
    level populations have not yet been established. During secondary emission, the absorption
    opacity \f$k^\text{abs}(\lambda)\f$ as a function of wavelength \f$\lambda\f$ is given by

    \f[ k^\text{abs}(\lambda) = \sum_\ell \, \frac{h c}{4\pi \lambda_\ell}(n_l\,B_{lu}-n_u\,B_{ul})
    \,\phi_\ell(\lambda) \f]

    where the sum with index \f$\ell\f$ runs over all supported lines, \f$u\f$ and \f$l\f$ are the
    indices of the energy levels before and after the corresponding transition, and
    \f$\phi_\ell(\lambda)\f$ is the Gaussian profile of line \f$\ell\f$.

    The calculation explicitly includes the Gaussian line profile caused by thermal motion of the
    molecules in the local bulk velocity frame of the cell. Moreover, the absorption profiles of
    all supported lines are superposed on top of each other. In practice, the implementation
    attempts to calculate just the terms that have a significant contribution at any given
    wavelength.

    */
class MolecularLineGasMix : public EmittingGasMix
{
    ITEM_CONCRETE(MolecularLineGasMix, EmittingGasMix,
                  "A gas mix supporting rotational and vibrational transitions in specific molecules and atoms")
        ATTRIBUTE_TYPE_INSERT(MolecularLineGasMix, "CustomMediumState,DynamicState")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[0")
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1000")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_DOUBLE_LIST(defaultCollisionPartnerRatios,
                             "the default relative abundances of the collisional partners")
        ATTRIBUTE_MIN_VALUE(defaultCollisionPartnerRatios, "[0")
        ATTRIBUTE_MAX_VALUE(defaultCollisionPartnerRatios, "1e50]")
        ATTRIBUTE_DEFAULT_VALUE(defaultCollisionPartnerRatios, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultCollisionPartnerRatios, "Level2")

        PROPERTY_DOUBLE(defaultMicroTurbulenceVelocity, "the default (non-thermal) micro-turbulence velocity")
        ATTRIBUTE_QUANTITY(defaultMicroTurbulenceVelocity, "velocity")
        ATTRIBUTE_MIN_VALUE(defaultMicroTurbulenceVelocity, "[0 km/s")
        ATTRIBUTE_MAX_VALUE(defaultMicroTurbulenceVelocity, "100000 km/s]")
        ATTRIBUTE_DEFAULT_VALUE(defaultMicroTurbulenceVelocity, "0 km/s")
        ATTRIBUTE_DISPLAYED_IF(defaultMicroTurbulenceVelocity, "Level2")

        PROPERTY_DOUBLE(maxChangeInLevelPopulations,
                        "the maximum relative change for the level populations in a cell to be considered converged")
        ATTRIBUTE_MIN_VALUE(maxChangeInLevelPopulations, "[0")
        ATTRIBUTE_MAX_VALUE(maxChangeInLevelPopulations, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxChangeInLevelPopulations, "0.05")
        ATTRIBUTE_DISPLAYED_IF(maxChangeInLevelPopulations, "Level2")

        PROPERTY_DOUBLE(maxFractionUnconvergedCells, "the maximum fraction of spatial cells that have not convergenced")
        ATTRIBUTE_MIN_VALUE(maxFractionUnconvergedCells, "[0")
        ATTRIBUTE_MAX_VALUE(maxFractionUnconvergedCells, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionUnconvergedCells, "0.001")
        ATTRIBUTE_DISPLAYED_IF(maxFractionUnconvergedCells, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the transition information for the configured specifies and determines
        and caches some further values used in the other functions. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns true, indicating that this material may have a negative absorption
        cross section. This happens when inverted level populations cause net stimulated emission.
        */
    bool hasNegativeExtinction() const override;

    /** This function returns true, indicating that the cross sections returned by this material
        mix depend on the values of specific state variables other than the number density. */
    bool hasExtraSpecificState() const override;

    /** This function returns DynamicStateType::PrimaryIfMergedIterations, indicating that this
        material mix has a dynamic medium state with updates that are considered to affect primary
        emission when the simulation has merged iterations, and only secondary emission if not. */
    DynamicStateType hasDynamicMediumState() const override;

    /** This function returns true, indicating that this material supports secondary line emission
        from gas. */
    bool hasLineEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this particular
        material mix as a list of SnapshotParameter objects. For this class, the function returns a
        descriptor for the number densities of the collisional partners and for the microturbulence
        level. Importing the kinetic gas temperature should be enabled through the corresponding
        standard configuration flag. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns a list
        containing descriptors for the properties defined in the input model and for a number of
        variables to hold the level populations derived from the radiation field and related
        information. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this material mix
        through the specificStateVariableInfo() function except for the number density of the
        species under consideration. For this class, the function uses the imported values, or if
        not available, the user-configured default values. The level populations are left at zero.
        */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** Based on the specified radiation field and the input model properties found in the given
        material state, this function determines the level populations for the supported
        transitions and stores these results back in the given material state. The function returns
        the update status as described for the UpdateStatus class. */
    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    /** This function returns true if the state of the medium component corresponding to this
        material mix can be considered to be converged based on the given spatial cell statistics,
        and false otherwise. The \em numCells, \em numUpdated and \em numNotConverged arguments
        specify respectively the number of spatial cells in the simulation, the number of cells
        updated during the latest update cycle, and the number of updated cells that have not yet
        converged. */
    bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a molecule or atom of the species under consideration. */
    double mass() const override;

    /** This function should return the absorption cross section using default properties. Because
        this value is hard to calculate for this material mix, this function returns zero. */
    double sectionAbs(double lambda) const override;

    /** This function should return the scattering cross section using default properties. Because
        this value is hard to calculate for this material mix, this function returns zero. */
    double sectionSca(double lambda) const override;

    /** This function should return the extinction cross section using default properties. Because
        this value is hard to calculate for this material mix, this function returns zero. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the absorption opacity \f$k^\text{abs}= n_\mathrm{mol}
        \varsigma^\text{abs}\f$ for the given wavelength and material state. The photon packet
        properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity \f$k^\text{sca}\f$ which is trivially zero at
        all wavelengths. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}\f$ for the given
        wavelength and material state, which equals the absorption opacity because the scattering
        opacity is zero. The photon packet properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the lines under consideration do not scatter. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the lines under consideration do not scatter. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

    /** This function returns a list with the line centers of the supported transitions. */
    Array lineEmissionCenters() const override;

    /** This function returns a list with the masses of the particles emitting each of the lines,
        i.e. the mass of a molecule or atom of the species under consideration in each case. */
    Array lineEmissionMasses() const override;

    /** This function returns a list with the line luminosities for the supported transitions in
        the spatial cell and medium component represented by the specified material state and the
        receiving material mix when it would be embedded in the specified radiation field. */
    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Temperature =======

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the temperature stored in the specific state for the relevant spatial
        cell and medium component. This value corresponds to the temperature defined by the input
        model at the start of the simulation. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // These data members are loaded from text files in setupSelfBefore()
    vector<double> _energyStates;  // the column of the energy states
    vector<double> _weights;       // the column of the weight of energy states
    vector<int> _indexUpRad;       // the column of the index of upper energy levels for the radiational transitions
    vector<int> _indexLowRad;      // the column of the index of lower energy levels for radiational transitions
    vector<double> _einsteinA;     // the column of Einstein A coefficients
    vector<int> _indexUpCol;       // the column of the index of upper energy levels for collisional transitions
    vector<int> _indexLowCol;      // the column of the index of lower energy levels for collisional transitions
    vector<vector<double>>
        _collisionCul;  // the column of collisional excitation coefficients for collisional transitions

    // These data members are precalculated in setupSelfBefore()
    vector<double> _einsteinBul;  // the column of Einstein Bul coefficients
    int _numLines{-1};
    vector<int> _indexRadTrans;  // index for the radiational transitions you use in this calculation
    vector<int> _indexColTrans;  // index for the collisional transitions you use in this calculation
    vector<int> _indexLines;     // index in the simulation's RF WLG for the bin containing rotational transition lines
};

////////////////////////////////////////////////////////////////////

#endif

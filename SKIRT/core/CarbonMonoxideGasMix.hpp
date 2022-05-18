/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CARBONMONOXIDEGASMIX_HPP
#define CARBONMONOXIDEGASMIX_HPP

#include "EmittingGasMix.hpp"

////////////////////////////////////////////////////////////////////

/** The CarbonMonoxideGasMix class describes the material properties related to selected rotational
    and rovibrational transitions in carbon monoxide (CO).

    The energy levels of CO are discretized into vibrational levels (quantum number \f$v\f$) and
    more subdivided rotational levels (quantum number \f$J\f$). The vibrational levels are
    quantized by the magnitude of the dipole moment, and the rotational levels are quantized by the
    magnitude of the dipole rotation. This class supports transitions involving vibrational levels
    \f$v=0\f$ and \f$v=1\f$ and rotational levels up to \f$J=20\f$. The selection rules allow two
    types of transitions between these energy levels. Rotational transitions change the rotational
    energy level by \f$\Delta J=\pm 1\f$. The corresponding lines are in the far-infrared to
    millimeter wavelength range, up to about 2600 \f$\mu\mathrm{m}\f$. Rovibrational transitions
    change both the rotational and vibrational levels at the same time (\f$\Delta J=\pm1\f$ and
    \f$v=0-1\f$). The corresponding lines, at least for the two supported vibrational levels, are
    all in the near infrared around 4.7 \f$\mu\mathrm{m}\f$.

    For each supported transition, the emission luminosity and self-absorption opacity in a given
    cell are determined from gas properties defined in the input model (CO number density,
    molecular hydrogen number density, kinetic gas temperature) and the local radiation field
    calculated by the simulation taking into account dust extinction. The class performs an
    iterative, non-LTE calculation including the effects of collisional transitions (excitation and
    de-excitation with molecular hydrogen as interaction partner) and photonic transitions
    (spontaneous emission, absorption, and induced emission). This allows establishing the energy
    level distribution (population) for a wide range of material densities and radiation fields.

    <b>Configuring the simulation</b>

    Simulations of the CO transition spectra should include primary sources and a dust medium in
    addition to a medium component configured with the carbon monoxide mix (this class).
    Furthermore, the \em simulationMode should be set to "DustAndGasEmission" and \em
    iterateSecondaryEmission should be enabled. The radiation field wavelength grid should properly
    resolve the UV, optical, and infrared wavelength range (relevant for dust emission) and the
    wavelength ranges of the supported CO transition lines.

    Separate instruments can be configured for the relevant wavelength ranges, e.g. using a
    logarithmic grid for the continuum spectrum and linear grids for the line profiles.

    During primary emission, the simulation determines the radiation field resulting from the
    primary sources and dust attenuation (CO line absorption is taken to be zero at this stage).
    The resulting radiation field allows a first estimation of the CO level populations for each
    spatial grid cell; this calculation happens in the updateSpecificState() function, which is
    invoked at the end of primary emission.

    During secondary emission, the simulation takes into account emission from all configured media
    (dust continuum and CO lines) in addition to absorption by these same media (dust and CO). For
    CO, the previously stored level populations allow calculating the line emission spectrum and
    the line self-absorption cross sections. This results in an updated radiation field, which will
    in turn influence the CO level populations (and for high optical depths, possibly the dust
    temperature), which in turn influences the secondary emission spectra. In order to obtain a
    self-consistent result, the simulation must therefore iterate over secondary emission.

    As indicated above, the input model must provide values for the spatial distribution of three
    medium properties: CO number density, molecular hydrogen number density, and kinetic gas
    temperature. These values remain constant during the simulation. Most often, this information
    will be read from an input file by associating the CO material mix with a subclass of
    ImportedMedium. For that medium component, the ski file attribute \em importTemperature
    <b>must</b> be set to 'true', and \em importMetallicity and \em importVariableMixParams must be
    left at 'false'. The additional column required by the CO material mix (molecular hydrogen
    number density) is automatically imported and is expected <b>after</b> all other columns. For
    example, if bulk velocities are also imported for this medium component (i.e. \em
    importVelocity is 'true'), the column order would be \f[ ..., n_\mathrm{CO}, T_\mathrm{kin},
    v_\mathrm{x}, v_\mathrm{y}, v_\mathrm{z}, n_\mathrm{H2} \f]

    For basic testing purposes, the CO material mix can also be associated with a geometric medium
    component. The geometry then defines the spatial density distribution (i.e.
    \f$n_\mathrm{CO}\f$), and the material mix offers configuration properties to specify a fixed
    default value for the gas temperature and for the number of hydrogen molecules per CO molecule
    that will be used across the spatial domain. In this case, the molecular hydrogen number
    density is thus implicitly defined based on the CO number density by a constant multiplier.

    <b>Level populations</b>

    The level populations \f$n_i\f$ (used later to calculate the emission luminosities and the
    absorption opacities) are obtained by solving the statistical equilibrium equations. The
    equation for the energy level with index \f$i\f$ (where indices increase from lowest to highest
    energy state) is given by

    \f[ \sum^{N_\mathrm{lp}}_{j>i} n_jA_{ji} + \sum^{N_\mathrm{lp}}_{j\neq i} \Big[ (n_jB_{ji} -
    n_iB_{ij})J_\lambda(\lambda_{ij})\Big] + \sum^{N_\mathrm{lp}}_{j\neq i}
    \big[n_jC_{ji}(n_\mathrm{H2} ,T_\mathrm{kin}) - n_iC_{ij}(n_\mathrm{H2},T_\mathrm{kin})\big]=0,
    \f]

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
    a thermal velocity corresponding to the kinetic temperature in the cell and the mass of a CO
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
class CarbonMonoxideGasMix : public EmittingGasMix
{
    ITEM_CONCRETE(CarbonMonoxideGasMix, EmittingGasMix,
                  "A gas mix supporting the carbon monoxide rotational and rovibrational transitions")
        ATTRIBUTE_TYPE_INSERT(CarbonMonoxideGasMix, "CustomMediumState,SemiDynamicState")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[3")  // gas temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_DOUBLE(defaultMolecularHydrogenRatio, "the default ratio of H2 over CO molecule numbers")
        ATTRIBUTE_MIN_VALUE(defaultMolecularHydrogenRatio, "[0")
        ATTRIBUTE_MAX_VALUE(defaultMolecularHydrogenRatio, "1e6]")
        ATTRIBUTE_DEFAULT_VALUE(defaultMolecularHydrogenRatio, "100")
        ATTRIBUTE_DISPLAYED_IF(defaultMolecularHydrogenRatio, "Level2")

        PROPERTY_DOUBLE(maxFractionNotConverged, "the maximum fraction of not-converged spatial cells")
        ATTRIBUTE_MIN_VALUE(maxFractionNotConverged, "[0")
        ATTRIBUTE_MAX_VALUE(maxFractionNotConverged, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionNotConverged, "0.001")
        ATTRIBUTE_DISPLAYED_IF(maxFractionNotConverged, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function determines and caches some values used in the other functions. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns true, indicating that the cross sections returned by this material
        mix depend on the values of specific state variables other than the number density. */
    bool hasExtraSpecificState() const override;

    /** This function returns true, indicating that this material has a semi-dynamic medium state.
        */
    bool hasSemiDynamicMediumState() const override;

    /** This function returns true, indicating that this material supports secondary line emission
        from gas. */
    bool hasLineEmission() const override;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this particular
        material mix as a list of SnapshotParameter objects. For this class, the function returns a
        descriptor for the molecular hydrogen density import parameter. Importing the kinetic gas
        temperature should be enabled through the corresponding standard configuration flag. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns a list
        containing descriptors for the properties defined in the input model (CO number density, H2
        number density, and temperature) and for a number of variables to hold the CO level
        populations derived from the radiation field when the semi-dynamic medium state is updated.
        */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this material mix
        through the specificStateVariableInfo() function except for the CO number density. For this
        class, the function initializes the H2 number density and the temperature to the specified
        imported values, or if not available, to the user-configured default values. The level
        populations are set to zero. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** Based on the specified radiation field and the input model properties found in the given
        material state, this function determines the level populations for the supported CO
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
    /** This function returns the mass of a neutral hydrogen atom. */
    double mass() const override;

    /** This function returns the absorption cross section per CO atom at the given wavelength and
        using the default gas properties configured for this material mix. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per CO atom, which is trivially zero for
        all wavelengths. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per CO atom at the given
        wavelength and using the default gas properties configured for this material mix. The
        extinction cross section is identical to the absorption cross section because the
        scattering cross section is zero. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the CO absorption opacity \f$k^\text{abs}= n_\mathrm{CO}
        \varsigma^\text{abs}\f$ for the given wavelength and material state. The photon packet
        properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the CO scattering opacity \f$k^\text{sca}\f$ which is trivially zero
        at all wavelengths. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the CO extinction opacity \f$k^\text{ext}=k^\text{abs}\f$ for the
        given wavelength and material state, which equals the absorption opacity because the
        scattering opacity is zero. The photon packet properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the CO lines do not scatter. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function does nothing because the CO lines do not scatter. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

    /** This function returns a list with the line centers of the supported CO transitions. */
    Array lineEmissionCenters() const override;

    /** This function returns a list with the masses of the particle emitting each of the lines,
        i.e. the CO molecule in each case. */
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
    // ...
};

////////////////////////////////////////////////////////////////////

#endif

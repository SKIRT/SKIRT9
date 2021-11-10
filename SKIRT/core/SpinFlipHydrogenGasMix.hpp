/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPINFLIPHYDROGENGASMIX_HPP
#define SPINFLIPHYDROGENGASMIX_HPP

#include "EmittingGasMix.hpp"

////////////////////////////////////////////////////////////////////

/** At some point in the future, the SpinFlipHydrogenGasMix class will hopefully describe the
    material properties related to the 21 cm spin-flip transition in neutral hydrogen, including
    emission and absorption. For now, it just serves as a stub for testing the framework that will
    enable this functionality.

    The spatial distributions for the gas density, metallicity, and temperature and for the local
    dust-to-gas ratio must be defined by the input model and are considered to be constant during
    the simulation. In this context, this material mix offers configuration properties to specify
    default values for these quantities that will be used by geometric media across the spatial
    domain. */
class SpinFlipHydrogenGasMix : public EmittingGasMix
{
    ITEM_CONCRETE(SpinFlipHydrogenGasMix, EmittingGasMix,
                  "A gas mix supporting the spin-flip 21 cm hydrogen transition")
        ATTRIBUTE_TYPE_INSERT(SpinFlipHydrogenGasMix, "CustomMediumState")

        PROPERTY_DOUBLE(defaultMetallicity, "the default metallicity of the gas")
        ATTRIBUTE_MIN_VALUE(defaultMetallicity, "]0")
        ATTRIBUTE_MAX_VALUE(defaultMetallicity, "0.2]")
        ATTRIBUTE_DEFAULT_VALUE(defaultMetallicity, "0.02")
        ATTRIBUTE_DISPLAYED_IF(defaultMetallicity, "Level2")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[3")  // gas temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e6]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_DOUBLE(defaultDustToGasRatio, "the default dust-to-gas ratio")
        ATTRIBUTE_MIN_VALUE(defaultDustToGasRatio, "]0")
        ATTRIBUTE_MAX_VALUE(defaultDustToGasRatio, "0.2]")
        ATTRIBUTE_DEFAULT_VALUE(defaultDustToGasRatio, "0.01")
        ATTRIBUTE_DISPLAYED_IF(defaultDustToGasRatio, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function determines the radiation field wavelength bin containing the UV field
        strength. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

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
        descriptor for the dust-to-gas ratio import parameter. Importing metallicity and
        temperature should be enabled through the corresponding standard configuration flags. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns a list
        containing descriptors for number density, metallicity, temperature, a custom variable to
        hold the dust-to-gas-ratio, and a custom variable to hold the UV field strength derived
        from the radiation field when the semi-dynamic medium state is updated. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this fragmented dust
        mix through the specificStateVariableInfo() function except for the number density. For
        this class, the function initializes the temperature, metallicity and dust-to-gas ratio to
        the specified imported values, or if not available, to the user-configured default values.
        The UV field strength is set to zero. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Medium state updates =======

    /** Based on the specified radiation field, the function obtains the UV field strength and
        stores it in the medium state for the specified cell. The function returns true if the
        medium state has indeed be changed, and false otherwise. */
    bool updateSpecificState(MaterialState* state, const Array& Jv) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a neutral hydrogen atom. */
    double mass() const override;

    /** This function returns the 21 cm absorption cross section per hydrogen atom at the given
        wavelength and using the default gas properties configured for this material mix. */
    double sectionAbs(double lambda) const override;

    /** This function returns the 21 scattering cross section per hydrogen atom, which is trivially
        zero for all wavelengths. */
    double sectionSca(double lambda) const override;

    /** This function returns the total 21 cm extinction cross section per hydrogen atom at the
        given wavelength and using the default gas properties configured for this material mix. The
        extinction cross section is identical to the absorption cross section because the
        scattering cross section is zero. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the 21 cm absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$
        for the given wavelength and material state. The photon packet properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the 21 cm scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$
        which is trivially zero at all wavelengths. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the 21 cm extinction opacity
        \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$ for the given wavelength and material state.
        The photon packet properties are not used. The extinction opacity equals the absorption
        opacity. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function throws a fatal error because the 21 cm line does not scatter. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w, Direction bfkobs,
                           Direction bfky, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function throws a fatal error because the 21 cm line does not scatter. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

    /** This function returns a list including a single item: the line center of the 21 cm hydrogen
        spinflip transition. */
    Array lineEmissionCenters() const override;

    /** This function returns a list including a single item: the mass of the particle emitting the
        21 cm line, i.e. the hydrogen atom. */
    Array lineEmissionMasses() const override;

    /** This function returns  a list including a single item: the 21 cm line luminosity
        in the spatial cell and medium component represented by the specified material state and
        the receiving material mix when it would be embedded in the specified radiation field. */
    Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const override;

    //======== Temperature =======

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the temperature stored in the specific state for the relevant spatial
        cell and medium component. Because the hydrogen temperature is not calculated
        self-consistently in our treatment, this value corresponds to the temperature defined by
        the input model at the start of the simulation. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    int _indexUV{-1};  // index in simulation's RF WLG for the bin containing 1000 A
};

////////////////////////////////////////////////////////////////////

#endif

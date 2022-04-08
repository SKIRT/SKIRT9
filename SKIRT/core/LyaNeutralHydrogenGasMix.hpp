/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYANEUTRALHYDROGENGASMIX_HPP
#define LYANEUTRALHYDROGENGASMIX_HPP

#include "DipolePhaseFunction.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** The LyaNeutralHydrogenGasMix class describes the material properties related to
    Lyman-alpha line transfer for a population of neutral hydrogen atoms, including support for
    polarization by scattering.

    The spatial distributions for both the mass density and the temperature of the neutral hydrogen
    gas must be defined by the input model and are considered to be constant during the simulation.
    In this context, this material mix offers a configuration property to specify a default gas
    temperature that is used by geometric media as a fixed temperature across the spatial domain.
    Imported media offer an \em importTemperature flag that allows defining a different gas
    temperature for each particle or cell. In this case, the default dust temperature configured
    for the material mix is ignored.

    This material mix also offers a configuration property to enable or disable support for
    polarization. It returns one of the \c Lya or \c LyaPolarization scattering modes depending on
    the configured value for this property. The dust-oriented functions for retrieving material
    properties are not used for the Lyman-alpha-specific aspects of the photon cycle, but some of
    them are used during setup, for example to normalize the mass of a geometric medium component
    based on optical depth. To this end, the scattering cross section is calculated using the
    configured default temperature and the absorption cross section is zero. */
class LyaNeutralHydrogenGasMix : public MaterialMix
{
    ITEM_CONCRETE(LyaNeutralHydrogenGasMix, MaterialMix, "neutral hydrogen for Lyman-alpha line transfer")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LyaNeutralHydrogenGasMix, "Lya")
        ATTRIBUTE_TYPE_INSERT(LyaNeutralHydrogenGasMix, "GasMix")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the neutral hydrogen gas")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[3")  // gas temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

        PROPERTY_BOOL(includePolarization, "include support for polarization")
        ATTRIBUTE_DEFAULT_VALUE(includePolarization, "false")
        ATTRIBUTE_DISPLAYED_IF(includePolarization, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function initializes the DipolePhaseFunction instance held by this class. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

    /** This function returns the value of the \em includePolarization flag, indicating whether the
        material mix supports polarization during scattering events or not. */
    bool hasPolarizedScattering() const override;

    /** This function returns true, indicating that scattering for the material mix is resonant. */
    bool hasResonantScattering() const override;

    /** This function returns true, indicating that the cross sections returned by this material
        mix depend on the values of specific state variables other than the number density, in this
        case the temperature. */
    bool hasExtraSpecificState() const override;

    /** This function returns true, indicating that a scattering interaction for this material mix
        may (and usually does) adjust the wavelength of the interacting photon packet. */
    bool hasScatteringDispersion() const override;

    //======== Medium state setup =======

public:
    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. See the description of the
        MaterialMix::specificStateVariableInfo() function for more information.

        The Lyman-alpha material mix requires a gas temperature in addition to the standard number
        density, so this function returns a list containing these two items. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes the specific state variables requested by this fragmented dust
        mix through the specificStateVariableInfo() function except for the number density. For the
        Lyman-alpha material mix, the function initializes the temperature to the specified
        imported temperature, or if this is not available, to the user-configured default
        temperature for this material mix. The metallicity and custom parameter arguments are
        ignored. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of neutral hydrogen atom. */
    double mass() const override;

    /** This function returns the Lyman-alpha absorption cross section per hydrogen atom
        \f$\varsigma^\text{abs}_\alpha(\lambda, T_\text{def})\f$ which is trivially zero for all
        wavelengths and temperatures. */
    double sectionAbs(double lambda) const override;

    /** This function returns the Lyman-alpha scattering cross section per hydrogen atom
        \f$\varsigma^\text{sca}_\alpha(\lambda, T_\text{def})\f$ at the given wavelength and using
        the default gas temperature configured for this material mix. */
    double sectionSca(double lambda) const override;

    /** This function returns the total Lyman-alpha extinction cross section per hydrogen atom
        \f$\varsigma^\text{ext}_\alpha(\lambda, T_\text{def})\f$ at the given wavelength and using
        the default gas temperature configured for this material mix. The extinction cross section
        is identical to the scattering cross section because the absorption cross section is zero.
        */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$,
        which is trivially zero for the Lyman-alpha material mix. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$ for
        the given wavelength and material state. The photon properties are not used. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$
        for the given wavelength and material state. The photon properties are not used. For the
        Lyman-alpha material mix, the extinction opacity equals the scattering opacity. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift,
        for the given wavelength, geometry, material state, and photon properties. See the
        description of the MaterialMix::peeloffScattering() function for more information.

        For the Lyman-alpha material mix, the function implements resonant scattering without or
        with support for polarization depending on the user-configured \em includePolarization
        property. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. For the Lyman-alpha material mix, the function implements resonant scattering
        without or with support for polarization depending on the user-configured \em
        includePolarization property. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Secondary emission =======

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the temperature stored in the specific state for the relevant spatial
        cell and medium component. Because the hydrogen temperature is not calculated
        self-consistently in our treatment, this value corresponds to the temperature defined by
        the input model at the start of the simulation. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // the dipole phase function helper instance - initialized during setup
    DipolePhaseFunction _dpf;
};

////////////////////////////////////////////////////////////////////

#endif

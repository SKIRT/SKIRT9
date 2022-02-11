/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONMIX_HPP
#define ELECTRONMIX_HPP

#include "ComptonPhaseFunction.hpp"
#include "DipolePhaseFunction.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ElectronMix class describes the material properties for a population of electrons.
    Electrons do not absorb photons. They do, however, significantly scatter photons. This process
    is described by Compton scattering, which converges to Thomson scattering at low photon
    energies. It is meaningful to implement both processes, because the calculations for Compton
    scattering are substantially slower than those for Thomson scattering.

    <b>Compton and Thomson scattering</b>

    For wavelengths shorter than 10 nm, this class models Compton scattering, which features a
    wavelength-dependent cross section and phase function, and which causes the photon energy
    (wavelength) to change during the interaction. This process is implemented through the
    ComptonPhaseFunction class; see that class for more information. The current implementation of
    Compton scattering does not support polarization.

    For wavelengths longer than 10nm, the scattering process can be described by elastic and
    wavelength-independent Thomson scattering. The scattering cross section is given by the
    well-known Thomson cross section (a constant) and the phase function is that of a dipole.
    Consequently, this class calls on the DipolePhaseFunction class to implement Thomson
    scattering; see that class for more information. In this wavelength regime, polarization is
    fully supported (if enabled by the user).

    The transition point between Compton and Thomson scattering can be justified as follows. For
    wavelengths much longer than 10 nm, the expression for the Compton cross section becomes
    numerically unstable. This could be solved by using a series expansion approximation at longer
    wavelengths. However, at a wavelength of 10 nm, the Compton cross section has approached the
    Thomson constant to within less than 0.05 per cent. It thus seems reasonable to transition to
    the much faster Thomson scattering at that wavelength point.

    <b>Thermal dispersion</b>

    By default, this class assumes that all electrons are at rest in the local frame (i.e., there
    is no movement other than the bulk velocity of the spatial cell containing the electrons). If
    the \em includeThermalDispersion configuration flag is enabled (and the simulation mode is
    panchromatic), a random thermal motion corresponding to the local temperature is added when
    performing a scattering interaction. The dispersion is not taken into account to determine the
    scattering cross section, because that value does not vary significantly for small wavelength
    shifts (and it does not vary at all for Thomson scattering).

    If the electron mix is associated with an ImportedMedium and the \em importTemperature flag for
    the medium is enabled, the local dispersion temperature is obtained from the imported file.
    Otherwise, the value configured for the \em defaultTemperature property is used instead,
    resulting in a constant temperature across space. */
class ElectronMix : public MaterialMix
{
    ITEM_CONCRETE(ElectronMix, MaterialMix, "a population of electrons")

        PROPERTY_BOOL(includePolarization, "include support for polarization")
        ATTRIBUTE_DEFAULT_VALUE(includePolarization, "false")
        ATTRIBUTE_DISPLAYED_IF(includePolarization, "Level2")

        PROPERTY_BOOL(includeThermalDispersion, "include thermal velocity dispersion")
        ATTRIBUTE_DEFAULT_VALUE(includeThermalDispersion, "false")
        ATTRIBUTE_RELEVANT_IF(includeThermalDispersion, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(includeThermalDispersion, "Level2")

        PROPERTY_DOUBLE(defaultTemperature, "the default temperature of the electron population")
        ATTRIBUTE_QUANTITY(defaultTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(defaultTemperature, "[3")    // temperature must be above local Universe T_CMB
        ATTRIBUTE_MAX_VALUE(defaultTemperature, "1e8]")  // higher temperatures cause relativistic dispersion
        ATTRIBUTE_DEFAULT_VALUE(defaultTemperature, "1e4")
        ATTRIBUTE_RELEVANT_IF(defaultTemperature, "Panchromatic&includeThermalDispersion")
        ATTRIBUTE_DISPLAYED_IF(defaultTemperature, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function initializes the DipolePhaseFunction instance held by this class. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Electrons. */
    MaterialType materialType() const override;

    /** This function returns the value of the \em includePolarization flag, indicating whether the
        material mix supports polarization during scattering events or not. */
    bool hasPolarizedScattering() const override;

    /** This function returns true when thermal dispersion is enabled and/or support for Compton
        scattering has been turned on (based on the simulation's wavelength range), indicating that
        in those cases a scattering interaction for the electron mix may (and usually does) adjust
        the wavelength of the interacting photon packet. If both thermal dispersion and support for
        Compton scattering are disabled, the function returns false. */
    bool hasScatteringDispersion() const override;

    //======== Medium state setup =======

public:
    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. See the description of the
        MaterialMix::specificStateVariableInfo() function for more information.

        For the electron mix class, the returned list always includes the specific state variable
        for number density and, in case the electron mix has thermal dispersion, it also includes
        the specific state variable for temperature. */
    vector<StateVariable> specificStateVariableInfo() const override;

    /** This function initializes any specific state variables requested by this material mix
        except for the number density. See the description of the
        MaterialMix::initializeSpecificState() function for more information.

        For this class, in case the electron mix has thermal dispersion, the function initializes
        the temperature to the specified imported temperature, or if this is not available, to the
        user-configured default temperature for this electron mix. */
    void initializeSpecificState(MaterialState* state, double metallicity, double temperature,
                                 const Array& params) const override;

    //======== Low-level material properties =======

public:
    /** This function returns the electron mass. */
    double mass() const override;

    /** This function returns the absorption cross section per electron
        \f$\varsigma^{\text{abs}}_{\lambda}\f$, which is trivially zero for all wavelengths
        \f$\lambda\f$. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per electron
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ which is constant and equal to the Thomson cross
        section for all wavelengths \f$\lambda\f$. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per electron
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ which is constant and equal to the Thomson cross
        section for all wavelengths \f$\lambda\f$. */
    double sectionExt(double lambda) const override;

    //======== High-level photon life cycle =======

    /** This function returns the absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$,
        which is trivially zero for electrons. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$ for
        the given material state. The wavelength and photon properties are not used, because the
        cross section is considered to be equal to the Thomson cross section for all wavelengths.
        */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$
        for the given material state. The wavelength and photon properties are not used, because
        the cross section is considered to be equal to the Thomson cross section for all
        wavelengths. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift
        for the given wavelength, geometry, material state, and photon properties. See the
        description of the MaterialMix::peeloffScattering() function for more information.

        For electrons, the function implements wavelenth-independent dipole scattering without or
        with support for polarization depending on the user-configured \em includePolarization
        property. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. For electrons, the function implements wavelenth-independent dipole
        scattering without or with support for polarization depending on the user-configured \em
        includePolarization property. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======================== Probing ========================

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. The implementation in this class ignores the radiation
        field and returns the temperature driving the thermal velocity dispersion for the medium
        component in the relevant spatial cell, or zero if thermal velocity dispersion is not
        enabled for this material mix. Because nothing in the simulation changes the electron
        temperature, this value corresponds to the temperature defined by the input model at the
        start of the simulation. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // flags initialized during setup
    bool _hasDispersion{false};  // true if thermal velocity dispersion is enabled
    bool _hasCompton{false};     // true if support for Compton scattering is enabled

    // the dipole and Compton phase function helper instances - initialized during setup
    DipolePhaseFunction _dpf;
    ComptonPhaseFunction _cpf;
};

////////////////////////////////////////////////////////////////////

#endif

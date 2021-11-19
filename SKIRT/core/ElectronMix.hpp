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

    Specifically, for wavelengths shorter than 10 nm, this class models Compton scattering, which
    features a wavelength-dependent cross section and phase function, and which causes the photon
    energy (wavelength) to change during the interaction. This process is implemented through the
    ComptonPhaseFunction class; see that class for more information. The current implementation of
    Compton scattering does not support polarization.

    For wavelengths longer than 10nm, the scattering process can be described by elastic and
    wavelength-independent Thomson scattering. The scattering cross section is given by the
    well-known Thomson cross section (a constant) and the phase function is that of a dipole.
    Consequently, this class calls on the DipolePhaseFunction class to implement Thomson
    scattering; see that class for more information. In this wavelength regime, polarization is
    fully supported (if enabled by the user).

    The transition point between Compton and Thomson scattering can be justified as follows. At a
    wavelength of 10 nm, the Compton cross section has approached the Thomson constant to within
    less than 0.05 per cent. At the same time, for wavelengths much longer than 10 nm, the
    expression for the Compton cross section becomes numerically unstable. Switching to Thomson
    scattering for longer wavelengths is thus important for accuracy in addition to performance.

    The current implementation of this class assumes that all electrons are at rest in the local
    bulk velocity frame (i.e they all move exactly at the local bulk velocity). */
class ElectronMix : public MaterialMix
{
    ITEM_CONCRETE(ElectronMix, MaterialMix, "a population of electrons")

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
        is MaterialType::Electrons. */
    MaterialType materialType() const override;

    /** This function returns the value of the \em includePolarization flag, indicating whether the
        material mix supports polarization during scattering events or not. */
    bool hasPolarizedScattering() const override;

    //======== Medium state setup =======

public:
    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. See the description of the
        MaterialMix::specificStateVariableInfo() function for more information.

        Electrons require just the standard specific state variable of type numberDensity , so this
        function returns a list containing a single item. */
    vector<StateVariable> specificStateVariableInfo() const override;

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
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, double w, Direction bfkobs,
                           Direction bfky, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. For electrons, the function implements wavelenth-independent dipole
        scattering without or with support for polarization depending on the user-configured \em
        includePolarization property. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======================== Data Members ========================

private:
    // flag indicating whether support for Compton scattering is enabled
    bool _hasCompton{false};

    // the dipole and Compton phase function helper instances - initialized during setup
    DipolePhaseFunction _dpf;
    ComptonPhaseFunction _cpf;
};

////////////////////////////////////////////////////////////////////

#endif

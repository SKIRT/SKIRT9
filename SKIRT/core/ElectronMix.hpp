/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONMIX_HPP
#define ELECTRONMIX_HPP

#include "DipolePhaseFunction.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ElectronMix class describes the material properties for a population of electrons,
    including support for polarization by scattering.

    Electrons do not absorb photons, and at the wavelengths relevant for SKIRT, scattering of
    photons by electrons can be described by elastic and wavelength-independent Thomson scattering.
    The scattering cross section is given by the well-known Thomson cross section (a constant) and
    the Mueller matrix elements and the phase function are those of a dipole. Consequently, this
    class calls on the DipolePhaseFunction class to implement most of its functionality. See there
    for more information. */
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
                           Direction bfky, const MaterialState* state, PhotonPacket* pp) const override;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. For electrons, the function implements wavelenth-independent dipole
        scattering without or with support for polarization depending on the user-configured \em
        includePolarization property. */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    /** This function returns an indicative temperature of the material mix when it would be
        embedded in a given radiation field. Because electrons don't absorb nor emit (within the
        range of physics supported here), the implementation in this class always returns zero. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

private:
    // the dipole phase function helper instance - initialized during setup
    DipolePhaseFunction _dpf;
};

////////////////////////////////////////////////////////////////////

#endif

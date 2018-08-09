/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONMIX_HPP
#define ELECTRONMIX_HPP

#include "MaterialMix.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The ElectronMix class describes the material properties for a population of electrons,
    including support for polarization by scattering.

    Electrons do not absorb photons, and at the wavelengths relevant for SKIRT, scattering of
    photons by electrons can be described by elastic and wavelength-independent Thomson scattering.
    The scattering cross section is given by the well-known Thomson cross section (a constant) and
    the Mueller matrix elements (and hence the phase function) can be expressed analytically as a
    function of just the scattering angle; see Bohren & Huffman (1998) or Wolf 2003 (Computer
    Physics Communications, 150, 99–115). */
class ElectronMix : public MaterialMix
{
    ITEM_CONCRETE(ElectronMix, MaterialMix, "a population of electrons")

    PROPERTY_BOOL(includePolarization, "include support for polarization")
        ATTRIBUTE_DEFAULT_VALUE(includePolarization, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches some pre-calculated arrays. */
    void setupSelfBefore() override;

    //======== Functionality levels =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Electrons. See the documentation of the MaterialMix class for more
        information. */
    MaterialType materialType() const override;

    /** This function returns the scattering mode supported by this material mix, which is
        ScatteringMode::MaterialPhaseFunction or ScatteringMode::SphericalPolarization depending on
        the value of the \em includePolarization flag. */
    ScatteringMode scatteringMode() const override;

    //======== Basic material properties =======

public:
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
        section for all wavelengths \f$\lambda\f$.. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ for the electron
        population, which is trivially equal to one for all wavelengths \f$\lambda\f$. */
    double albedo(double lambda) const override;

    //======== Scattering with material phase function =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angle cosine \f$\cos\theta\f$, where the phase function is normalized as \f[\int_{-1}^1
        \Phi_\lambda(\cos\theta) \,\mathrm{d}\cos\theta =2.\f] */
    double phaseFunctionValueForCosine(double lambda, double costheta) const override;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$. */
    double generateCosineFromPhaseFunction(double lambda) const override;

    //======== Polarization through scattering by spherical particles =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angles \f$\theta\f$ and \f$\phi\f$, and for the specified incoming polarization state. The
        phase function is normalized as \f[\int\Phi_\lambda(\theta,\phi) \,\mathrm{d}\Omega
        =4\pi.\f] */
   double phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const override;

    /** This function generates random scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from
        the phase function \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$, and for the
        specified incoming polarization state. The results are returned as a pair of numbers in the
        order \f$\theta\f$ and \f$\phi\f$. */
    std::pair<double,double> generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const override;

    /** This function applies the Mueller matrix transformation for the specified wavelength
        \f$\lambda\f$ and scattering angle \f$\theta\f$ to the given polarization state (which
        serves as both input and output for the function). */
    void applyMueller(double lambda, double theta, StokesVector* sv) const override;

    //======================== Data Members ========================

private:
    // precalculated discretizations - initialized during setup
    Array _phiv;                    // indexed on f
    Array _phi1v;                   // indexed on f
    Array _phisv;                   // indexed on f
    Array _phicv;                   // indexed on f
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COMPTONPHASEFUNCTION
#define COMPTONPHASEFUNCTION

#include "Array.hpp"
#include "Direction.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** The ComptonPhaseFunction helper class represents Compton scattering of photons by free
    electrons. The current implementation does not support polarization.

    Compton scattering forms the extension of Thomson scattering into the high-photon-energy
    domain. The scattering event is elastic (i.e. the photon wavelength is changed by the
    interaction) and its cross section and phase function depend on the wavelength of the incoming
    photon. Compared to Thomson scattering, the total cross section is reduced at higher energies.
    With \f$x = (h\nu)/(m_e c^2)\f$ the incoming photon energy scaled to the electron rest energy,
    the ratio of the total Compton scattering cross section \f$\sigma_\mathrm{C}(x)\f$ over the
    (constant) total Thomson cross section \f$\sigma_\mathrm{T}\f$ is given by:

    \f[ \frac{\sigma_\mathrm{C}(x)}{\sigma_\mathrm{T}} = \frac{3}{4} \bigg( \frac{1+x}{(1+2x)^2} +
    \frac{2}{x^2} + \Big(\frac{1}{2x}-\frac{1+x}{x^3}\Big) \cdot \ln(1+2x) \bigg). \f]

    For an interaction with incoming photon energy \f$x\f$ and scattering angle \f$\theta\f$, the
    energy change \f$E_\mathrm{out}/E_\mathrm{in}\f$ is given by the Compton factor \f$C(x, \theta)
    = \big(1+x(1-\cos \theta)\big)^{-1}\f$. Equivalently, the wavelength change is given by
    \f$\lambda_\mathrm{out} / \lambda_\mathrm{in} = \big(1+x(1-\cos \theta)\big)\f$.

    */
class ComptonPhaseFunction
{
    //============= Construction - Setup - Destruction =============

public:
    /** This function caches a pointer to the simulation's random number generator and
        pre-calculates some data used for sampling from the relevant phase function. The function
        must be called during setup (i.e. in single-thread mode). */
    void initialize(Random* random);

    //======== Absorption =======

    /** This function calculates and returns the Compton scattering cross section for the given
        wavelength, normalized to the (constant) Thomson cross section. */
    double sectionSca(double lambda) const;

    //======== Scattering =======

private:
    /** This function returns the value of the scattering phase function \f$\Phi(x, \cos\theta)\f$
        for the specified incoming photon energy \f$x\f$ and scattering angle cosine
        \f$\cos\theta\f$, where the phase function is normalized for all \f$x\f$ as \f[\int_{-1}^1
        \Phi(x, \cos\theta) \,\mathrm{d}\cos\theta =2.\f] */
    double phaseFunctionValueForCosine(double x, double costheta) const;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi(x, \cos\theta)\f$ for the specified incoming photon energy \f$x\f$. */
    double generateCosineFromPhaseFunction(double x) const;

public:
    /** This function calculates the contribution of a Compton scattering event to the peel-off
        photon luminosity for the given geometry and wavelength, and determines the adjusted
        wavelength of the outgoing photon packet. The luminosity contribution is added to the
        incoming value of the \em I argument, and the adjusted wavelength is stored in the \em
        lambda argument. The \em w argument specifies the relative opacity weighting factor for
        this medium component. See the description of the MaterialMix::peeloffScattering() function
        for more information. */
    void peeloffScattering(double& I, double& lambda, double w, Direction bfk, Direction bfkobs) const;

    /** Given the incoming photon packet wavelength and direction this function calculates a
        randomly sampled new propagation direction for a Compton scattering event, and determines
        the adjusted wavelength of the outgoing photon packet. The adjusted wavelength is stored in
        the \em lambda argument, and the direction is returned. */
    Direction performScattering(double& lambda, Direction bfk) const;

    //======================== Data Members ========================

private:
    // the simulation's random number generator - initialized during construction
    Random* _random{nullptr};

    // precalculated discretizations - initialized during construction
};

////////////////////////////////////////////////////////////////////

#endif

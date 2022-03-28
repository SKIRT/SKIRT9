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

    The normalized phase function for an interaction with incoming photon energy \f$x\f$ is given
    by:

    \f[ \Phi(x, \theta) = \frac{\sigma_\mathrm{T}}{\sigma_\mathrm{C}(x)} \, \frac{3}{4}
    \left[ C^3(x, \theta) + C(x, \theta) -C^2(x, \theta)\sin^2\theta \right], \f]

    where \f$\theta\f$ is the scattering angle and \f$C(x, \theta)\f$ is the Compton factor defined
    earlier.

    <b>Sampling from the phase function</b>

    To draw a random scattering angle from the phase function, we use the algorithm described by
    Hua et al. 1997 (Computers in Physics 11, 660), which is a variation of the technique first
    suggested by Pei 1979 and often referred to as Khan's technique. A combination of composition
    and rejection methods, the algorithm avoids expensive operations and has a rejection rate of
    about 1/3 depending on the energy.

    Using our notation for the scaled energy \f$x\f$ of the incoming photon (see above), Hua et al.
    1997 define the doubled scaled incoming photon energy \f$\epsilon=2 x\f$ and the inverse
    Compton factor \f$r = 1 + x (1-\cos\theta)\f$. The sampling algorithm draws a random number for
    \f$r\f$, i.e. from the probability distribution for the inverse Compton factor at a given
    energy. The scattering angle can then easily be obtained from the definition of the inverse
    Compton factor. */
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
        lambda argument. See the description of the MaterialMix::peeloffScattering() function
        for more information. */
    void peeloffScattering(double& I, double& lambda, Direction bfk, Direction bfkobs) const;

    /** Given the incoming photon packet wavelength and direction this function calculates a
        randomly sampled new propagation direction for a Compton scattering event, and determines
        the adjusted wavelength of the outgoing photon packet. The adjusted wavelength is stored in
        the \em lambda argument, and the direction is returned. */
    Direction performScattering(double& lambda, Direction bfk) const;

    //======================== Data Members ========================

private:
    // the simulation's random number generator - initialized by initialize()
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif

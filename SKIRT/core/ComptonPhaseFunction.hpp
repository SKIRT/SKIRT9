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

/** The ComptonPhaseFunction helper class represents the Compton scattering process. The current
    implementation does not support polarization. */
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
    /** This function returns the value of the scattering phase function \f$\Phi(\cos\theta)\f$ for
        the specified scattering angle cosine \f$\cos\theta\f$, where the phase function is
        normalized as \f[\int_{-1}^1 \Phi(\cos\theta) \,\mathrm{d}\cos\theta =2.\f] */
    double phaseFunctionValueForCosine(double costheta) const;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi(\cos\theta)\f$. */
    double generateCosineFromPhaseFunction() const;

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

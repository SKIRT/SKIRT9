/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYUTILS_HPP
#define LYUTILS_HPP

#include "Direction.hpp"

class Configuration;
class Random;

////////////////////////////////////////////////////////////////////

/** This namespace offers utility functions for the resonant scattering of photons by a specific
    atomic or ionic transition, such as the hydrogen Lyman-alpha line or one of the many X-ray
    resonance lines supported by the XRayIonicGasMix class. The functions are generic in the
    transition's rest-frame wavelength, Einstein A-coefficient, and statistical weight, which are
    passed in as arguments rather than hard-coded.

    When a photon resonant with a transition is absorbed by an ion in the transition's lower
    level, the ion is excited to the upper level and a new photon of (nearly) the same rest-frame
    wavelength is emitted almost immediately as a result of the subsequent downward transition.
    This happens fast enough that we can consider the combined process as a single scattering
    event.

    <b>Scattering cross section</b>

    The cross section for resonant scattering of a photon by a single ion can be derived using
    quantum mechanical considerations, resulting in a sharply peaked profile as a function of the
    photon wavelength in the ion's rest frame. Because each ion has its own velocity, a photon with
    a given wavelength in the gas rest frame will appear Doppler shifted to a slightly different
    wavelength for each ion in the gas. To compute the cross section for a collection of moving
    ions, we must therefore convolve the single-ion cross section with the ion velocity
    distribution, which in turn depends on the gas temperature.

    Assuming a Maxwell-Boltzmann velocity distribution, we define the thermal Doppler width
    \f$v_\mathrm{th}\f$ as \f[ v_\mathrm{th} = \sqrt{\frac{2 k_\mathrm{B} T}{m}} \f] where
    \f$k_\mathrm{B}\f$ is the Boltzmann constant, \f$m\f$ is the mass of the scattering ion, and
    \f$T\f$ is the temperature of the gas. We then introduce the dimensionless frequency variable
    \f$x\f$, defined as \f[ x = \frac{\nu - \nu_0}{\nu_0} \,\frac{c}{v_\mathrm{th}} =
    \frac{v_\mathrm{p}}{v_\mathrm{th}} \f] where \f$\nu=c/\lambda\f$ is the regular frequency
    variable, \f$\nu_0=c/\lambda_0\f$ is the frequency at the transition's rest-frame center,
    \f$\lambda_0\f$ is the corresponding rest-frame wavelength, and \f$c\f$ is the speed of light
    in vacuum. The last equality introduces the velocity shift \f$v_\mathrm{p}\f$ of the photon
    frequency relative to the line center, defined by \f[ \frac{v_\mathrm{p}}{c} = \frac{\nu -
    \nu_0}{\nu_0} \approx -\frac{\lambda - \lambda_0}{\lambda_0} \f] where the approximate equality
    holds for \f$v_\mathrm{p}\ll c\f$.

    After neglecting some higher order terms, the convolution of the single-ion profile with the
    Maxwell-Boltzmann velocity distribution yields the following expression for the
    velocity-weighted scattering cross section \f$\sigma(x)\f$ as a function of the dimensionless
    photon frequency \f$x\f$: \f[ \sigma(x) = \sigma_0\,H(a,x) \f] where the cross section at the
    line center \f$\sigma_0\f$ is given by \f[ \sigma_0 = \frac{g\,\lambda_0^3\,A}{8\,\pi^{3/2}\,
    v_\mathrm{th}}, \f] with \f$g\f$ the statistical weight of the upper level and \f$A\f$ the
    Einstein A-coefficient of the transition; the Voigt parameter \f$a\f$ is given by \f[ a =
    \frac{\lambda_0\,\Gamma}{4\pi\,v_\mathrm{th}} \f] with \f$\Gamma\f$ the natural line width; and
    the Voigt function \f$H(a,x)\f$ is defined and evaluated by the VoigtProfile class, which also
    provides the sampling algorithm used below.

    <b>Frequency shift due to atom velocity</b>

    In most astrophysical conditions, the energy of the resonant photon before and after
    scattering is identical in the frame of the interacting ion. This is because the life time of
    the upper level is very short so that it is not perturbed over this short time interval.
    Because of the random thermal motion of the ion, energy conservation in the ion's frame
    translates to a change in the energy of the incoming and outgoing photon that depends on the
    velocity of the ion and the scattering direction. Given the velocity of the ion \f$\bf{v}\f$,
    we define the dimensionless ion velocity as \f${\bf{u}}={\bf{v}}/v_\mathrm{th}\f$. Denoting the
    propagation direction and dimensionless frequency of the photon before (after) scattering with
    \f$\bf{k}_\mathrm{in}\f$ and \f$x_\mathrm{in}\f$ (\f$\bf{k}_\mathrm{out}\f$ and
    \f$x_\mathrm{out}\f$), the resulting frequency change can be written as \f[x_\mathrm{out} =
    x_\mathrm{in} - {\bf{u}}\cdot{\bf{k}}_\mathrm{in} + {\bf{u}}\cdot{\bf{k}}_\mathrm{out} \f] This
    analysis ignores the energy transferred from the photon to the ion through recoil, an
    approximation that is justified in regular astrophysical conditions.

    Assuming a Maxwell-Boltzmann velocity distribution for the ions, the two components of the
    dimensionless ion velocity \f$\bf{u}\f$ that are orthogonal to the incoming photon direction
    \f$\bf{k}_\mathrm{in}\f$ have a Gaussian probability distribution with zero mean and a standard
    deviation of \f$1/\sqrt{2}\f$. The parallel component is more complicated: its probability
    distribution given \f$x_\mathrm{in}\f$ is proportional to both the Gaussian ion velocity
    distribution and the single-ion scattering cross section, reflecting the preference for
    photons to be scattered by ions to which they appear close to resonance; sampling this
    component is delegated to VoigtProfile::sample().

    Unlike a treatment specific to a single, fixed transition, the functions in this namespace do
    not themselves select between the isotropic and dipole scattering phase functions -- since
    that choice depends on the angular-momentum quantum numbers of the transition being scattered,
    which vary between the resonance lines supported by a caller such as XRayIonicGasMix, the
    caller determines the phase function on its own (for example from the transition's total
    angular momentum) and uses sampleAtomVelocity() below only for the atom velocity itself. */
namespace LyUtils
{
    /** This function returns the resonant scattering cross section \f$\sigma(\lambda)\f$ for a
        single transition, given the photon wavelength \em lambda in the local gas rest frame, the
        transition's rest-frame central wavelength \em center, Einstein A-coefficient \em A, and
        upper-level statistical weight \em g, and the Doppler width \em vth and Voigt parameter
        \em a corresponding to the gas temperature and the mass of the scattering ion, all as
        defined in the documentation above. The \em vth argument must be the Doppler width
        \f$\sqrt{2k_\mathrm{B}T/m}\f$, not the thermal velocity \f$\sqrt{k_\mathrm{B}T/m}\f$, and
        \em a must be constructed using that same \em vth (see the documentation above). */
    double section(double lambda, double center, double vth, double A, double a, double g);

    /** This function draws a random ion velocity as seen by an incoming photon resonant with a
        transition with rest-frame central wavelength \em center, Doppler width \em vth, and Voigt
        parameter \em a (all as defined in the documentation above), from the appropriate
        probability distributions, reflecting the preference for photons to be scattered by ions to
        which they appear close to resonance. Unlike a treatment specific to a single transition,
        this function does not select between the isotropic and dipole phase functions -- that
        decision, if needed, is left to the caller (see the documentation above).

        The \em lambda argument specifies the photon packet wavelength as it is perceived in the
        local gas frame. The \em T and \em nH arguments specify the gas temperature and the number
        density of the scattering species in the current spatial cell; together with the globally
        configured acceleration scheme (exposed as Configuration::lyaAccelerationScheme() and
        Configuration::lyaAccelerationStrength(), named for its original Lyman-alpha context but
        applied here regardless of which transition is being sampled), these optionally bias the
        sampled velocity towards the line wings so as to reduce the number of scattering events
        needed to escape a very optically thick medium. Note that \em T and \em vth both relate to
        the same physical temperature but are independent arguments because they enter the
        calculation differently: \em vth (specific to the scattering ion's mass) sets the width of
        the dimensionless frequency scale, while \em T enters the acceleration scheme's own scaling
        relation directly.

        The function proceeds as follows:

        - Convert the perceived photon packet wavelength to the corresponding dimensionless
        frequency using its definition given above.

        - Draw values for the components of the dimensionless ion velocity parallel and orthogonal
        to the incoming photon packet from the probability distribution described above (delegating
        to VoigtProfile::sample() for the parallel component) and from Gaussian distributions,
        respectively, optionally biased according to the configured acceleration scheme.

        - Transform the dimensionless frequency into the rest frame of the ion as described above.

        - Multiply the dimensionless ion velocity by the Doppler width to obtain the actual
        physical ion velocity, and return it. */
    Vec sampleAtomVelocity(double lambda, double center, double vth, double a, double T, double nH, Direction kin,
                           Configuration* config, Random* random);

    /** This function returns the Doppler-shifted wavelength in the gas bulk rest frame after a
        resonant scattering event, given the incoming wavelength in the gas bulk rest frame, the
        velocity of the interacting ion, and the incoming and outgoing photon packet directions.
        */
    double shiftWavelength(double lambda, const Vec& vatom, const Direction& kin, const Direction& kout);
}

////////////////////////////////////////////////////////////////////

#endif

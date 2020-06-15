/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYAUTILS_HPP
#define LYAUTILS_HPP

#include "Direction.hpp"
#include "Range.hpp"
class Configuration;
class Random;

////////////////////////////////////////////////////////////////////

/** This namespace offers utility functions related to the treatment of Lyman-alpha (Lya) line
    transfer.

    When a Lya photon is absorbed by a neutral hydrogen atom in the ground state, the atom is
    excited to the 2p energy level and a new Lya photon is emitted almost immediately as a result
    of the subsequent downward transition. This happens fast enough that we can consider the
    combined process as a scattering event.

    <b>Scattering cross section</b>

    The cross section for a single hydrogen atom of the Lya scattering of a photon can be derived
    using quantum mechanical considerations, resulting in a sharply peaked profile as a function of
    the photon wavelength in the atom's rest frame. Because each atom has its own velocity, a
    photon with a given wavelength in the gas rest frame will appear Doppler shifted to a slightly
    different wavelength for each atom in the gas. To compute the Lya absorption cross section for
    a collection of moving atoms, we must therefore convolve the single-atom cross section with the
    atom velocity distribution, which in turn depends on the gas temperature.

    Assuming a Maxwell-Boltzmann velocity distribution, we define the characteristic thermal
    velocity \f$v_\mathrm{th}\f$ as \f[ v_\mathrm{th} = \sqrt{\frac{2 k_\mathrm{B}
    T}{m_\mathrm{p}}} \f] where \f$k_\mathrm{B}\f$ is the Boltzmann constant, \f$m_\mathrm{p}\f$ is
    the proton mass, and \f$T\f$ is the temperature of the gas. We then introduce the dimensionless
    frequency variable \f$x\f$, defined as \f[ x = \frac{\nu - \nu_\alpha}{\nu_\alpha}
    \,\frac{c}{v_\mathrm{th}} = \frac{v_\mathrm{p}}{v_\mathrm{th}} \f] where \f$\nu=c/\lambda\f$ is
    the regular frequency variable, \f$\nu_\alpha=c/\lambda_\alpha\f$ is the frequency at the Lya
    line center, \f$\lambda_\alpha\f$ is the wavelength at the Lya line center, and \f$c\f$ is the
    speed of light in vacuum. The last equality introduces the velocity shift \f$v_\mathrm{p}\f$ of
    the photon frequency relative to the Lya line center, defined by \f[ \frac{v_\mathrm{p}}{c} =
    \frac{\nu - \nu_\alpha}{\nu_\alpha} \approx -\frac{\lambda - \lambda_\alpha}{\lambda_\alpha}
    \f] where the approximate equality holds for \f$v_\mathrm{p}\ll c\f$.

    After neglecting some higher order terms, the convolution of the single-atom profile with the
    Maxwell-Boltzmann velocity distribution yields the following expression for the
    velocity-weighted Lya scattering cross section \f$\sigma_\alpha(x,T)\f$ of a hydrogen atom in
    gas at temperature \f$T\f$ as a function of the dimensionless photon frequency \f$x\f$: \f[
    \sigma_\alpha(x,T) = \sigma_{\alpha,0}(T)\,H(a_\mathrm{v}(T),x) \f] where the cross section at
    the line center \f$\sigma_{\alpha,0}(T)\f$ is given by \f[ \sigma_{\alpha,0}(T) =
    \frac{3\lambda_\alpha^2 a_\mathrm{v}(T)}{2\sqrt{\pi}}, \f] the Voigt parameter
    \f$a_\mathrm{v}(T)\f$ is given by \f[ a_\mathrm{v}(T) = \frac{A_\alpha}{4\pi\nu_\alpha}
    \,\frac{c}{v_\mathrm{th}} \f] with \f$A_\alpha\f$ the Einstein A-coefficient of the Lya
    transition; and the Voigt function \f$H(a_\mathrm{v},x)\f$ is defined by \f[ H(a_\mathrm{v},x)
    = \frac{a_\mathrm{v}}{\pi} \int_{-\infty}^\infty \frac{\mathrm{e}^{-y^2} \,\mathrm{d}y}
    {(y-x)^2+a_\mathrm{v}^2} \f] which is normalized so that \f$H(a_\mathrm{v},0) = 1\f$.

    <b>Frequency shift due to atom velocity</b>

    In most astrophysical conditions, the energy of the Lyman-alpha photon before and after
    scattering is identical in the frame of the interacting atom. This is because the life time of
    the atom in its 2p state is very short so that it is not perturbed over this short time
    interval. Because of the random thermal motion of the atom, energy conservation in the atom's
    frame translates to a change in the energy of the incoming and outgoing photon that depends on
    the velocity of the atom and the scattering direction. Given the velocity of the atom
    \f$\bf{v}\f$, we define the dimensionless atom velocity as
    \f${\bf{u}}={\bf{v}}/v_\mathrm{th}\f$. Denoting the propagation direction and dimensionless
    frequency of the photon before (after) scattering with \f$\bf{k}_\mathrm{in}\f$ and
    \f$x_\mathrm{in}\f$ (\f$\bf{k}_\mathrm{out}\f$ and \f$x_\mathrm{out}\f$), the resulting
    frequency change can be written as \f[x_\mathrm{out} = x_\mathrm{in} -
    {\bf{u}}\cdot{\bf{k}}_\mathrm{in} + {\bf{u}}\cdot{\bf{k}}_\mathrm{out} \f] This analysis
    ignores the energy transferred from the photon to the atom through recoil, an approximation
    that is justified in regular astrophysical conditions.

    Assuming a Maxwell-Boltzmann velocity distribution for the atoms, the two components of the
    dimensionless atom velocity \f$\bf{u}\f$ that are orthogonal to the incoming photon direction
    \f$\bf{k}_\mathrm{in}\f$ have a Gaussian probability distribution with zero mean and a standard
    deviation of \f$1/\sqrt{2}\f$. The parallel component is more complicated. We denote the
    dimensionless atom velocity component parallel to the incoming photon direction as
    \f$u_\parallel\f$. The probability distribution \f$P(u_\parallel|x_\mathrm{in})\f$ for this
    component is proportional to both the Gaussian atom velocity distribution and the Lyman-alpha
    scattering cross section for a single atom, reflecting the preference for photons to be
    scattered by atoms to which they appear close to resonance. With the proper normalization this
    leads to \f[ P(u_\parallel|x_\mathrm{in}) = \frac{a_\mathrm{v}}{\pi
    H(a_\mathrm{v},x_\mathrm{in})}\,
    \frac{\mathrm{e}^{-u_\parallel^2}}{(u_\parallel-x_\mathrm{in})^2+a_\mathrm{v}^2} \f]

    <b>Scattering phase function</b>

    Lyman-alpha scattering takes one of two forms: isotropic scattering or dipole scattering; the
    latter is also called Rayleigh scattering. The corresponding phase functions depend only on the
    cosine of the scattering angle \f$\mu={\bf{k}}_\mathrm{in}\cdot{\bf{k}}_\mathrm{out}\f$. With
    normalization \f$\int_{-1}^1 P(\mu)\,\mathrm{d}\mu = 2\f$ these phase functions can be written,
    respectively, as \f[ P(\mu) = 1 \f] and \f[ P(\mu) = \frac{3}{4}(\mu^2+1) \f]

    Quantum mechanical considerations lead to a simple recipe for selecting the appropriate phase
    function depending on whether the incoming photon frequency is in the core or in the wings of
    the cross section. The recipe prescribes to treat 1/3 of all core scattering events as dipole,
    and the remaining 2/3 as isotropic; and to treat all wing scattering events as dipole. For the
    purpose of this recipe, the scattering event is considered to occur in the core if the incoming
    dimensionless photon frequency in the rest frame of the interacting atom is smaller than a
    critical value, \f$|x|<0.2\f$. */
namespace LyaUtils
{
    /** This function returns the Lyman-alpha scattering cross section per hydrogen atom
        \f$\sigma_\alpha(\lambda, T)\f$ at the given photon wavelength and gas temperature, using
        the definition given in the class header. */
    double section(double lambda, double T);

    /** This function draws a random hydrogen atom velocity as seen by an incoming photon from the
        appropriate probability distributions, reflecting the preference for photons to be
        scattered by atoms to which they appear close to resonance. In addition, it determines
        whether the photon scatters through the isotropic or dipole phase function.

        The function arguments include the photon packet wavelength as it is perceived in the local
        gas frame and the hydrogen temperature and number density in the current spatial cell. The
        latter two values are used in the variable acceleration scheme. The return value is a pair:
        the first item is the atom velocity and the second item is true for the dipole phase
        function and false for isotropic scattering.

        The function proceeds as follows:

        - Convert the perceived photon packet wavelength to the corresponding dimensionless
        frequency using its definition given in the class header.

        - Draw values for the components of the dimensionless atom velocity parallel and orthogonal
        to the incoming photon packet from the probability distribution described in the class
        header and from Gaussian distributions, respectively.

        - Transform the dimensionless frequency into the rest frame of the atom as described in the
        class header.

        - Select the isotropic or the dipole phase functions according to the scheme discussed in
        the class header, using the dimensionless frequency in the atom frame.

        - Multiply the dimensionless atom velocity by the thermal velocity to obtain the actual
        physical atom velocity.

        - Return the atom velocity and a flag indicating the selected phase function.

        */
    std::pair<Vec, bool> sampleAtomVelocity(double lambda, double T, double nH, Direction kin, Configuration* config,
                                            Random* random);

    /** This function returns the Doppler-shifted wavelength in the gas bulk rest frame after a
        Lyman-alpha scattering event, given the incoming wavelength in the gas bulk rest frame, the
        velocity of the interacting atom, and the incoming and outgoing photon packet directions.
        */
    double shiftWavelength(double lambda, const Vec& vatom, const Direction& kin, const Direction& kout);
}

////////////////////////////////////////////////////////////////////

#endif

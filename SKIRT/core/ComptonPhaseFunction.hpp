/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COMPTONPHASEFUNCTION
#define COMPTONPHASEFUNCTION

#include "Array.hpp"
#include "Direction.hpp"
class Random;
class StokesVector;

////////////////////////////////////////////////////////////////////

/** The ComptonPhaseFunction helper class represents Compton scattering of photons by free
    electrons, with optional support for polarization by scattering.

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

    <b>Unpolarized Compton scattering</b>

    The normalized phase function for an interaction with incoming photon energy \f$x\f$ is given
    by:

    \f[ \Phi(x, \theta) = \frac{\sigma_\mathrm{T}}{\sigma_\mathrm{C}(x)} \, \frac{3}{4}
    \left[ C^3(x, \theta) + C(x, \theta) -C^2(x, \theta)\sin^2\theta \right], \f]

    where \f$\theta\f$ is the scattering angle and \f$C(x, \theta)\f$ is the Compton factor defined
    earlier.

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
    Compton factor.

    <b>Polarized Compton scattering</b>

    The Müller matrix describing Compton scattering can be expressed as a function of the
    scattering angle \f$\theta\f$ and the incoming photon energy \f$x\f$ as follows (Fano 1949):

    \f[ {\bf{M}}(x,\theta) \propto \begin{pmatrix} C^3(x,\theta) + C(x,\theta) -C^2(\theta,
    x)\sin^2\theta & -C^2(x,\theta)\sin^2\theta & 0 & 0 \\ -C^2(x,\theta)\sin^2\theta &
    C^2(x,\theta) (1+\cos^2\theta) & 0 & 0 \\ 0 & 0 & 2C^2(x,\theta)\cos\theta & 0 \\ 0 & 0 & 0 &
    (C^3(x,\theta) + C(x,\theta)) \cos\theta \end{pmatrix}, \f]

    with \f$C(x, \theta)\f$ the Compton factor defined earlier, and assuming a random distribution
    for the electron spin direction, which causes the non-diagonal terms in the forth row and forth
    column to be zero (Depaola 2003). This matrix has five independent coefficients and converges
    to the Müller matrix for Thomson scattering at low photon energies, i.e. for \f$x\to 0\f$ and
    thus \f$C(x,\theta) \to 1\f$.

    Following the procedure of Peest et al. 2017 Sect 3.3 starting from the above Müller matrix,
    we obtain the following expression for the normalized phase function:

    \f[ \Phi(x, \theta, \varphi, {\bf{S}}) = \frac{\sigma_\mathrm{T}}{\sigma_\mathrm{C}(x)} \,
    \frac{3}{4} \left[{\rm S}_{11}(x,\theta) + {\rm S}_{12}(x,\theta)
    \,P_\text{L}\cos2(\varphi-\gamma)\right], \f]

    with \f$\theta\f$ and \f$\varphi\f$ the inclination and azimuth angle of the scattering
    geometry, \f$\sigma_\mathrm{C}(x)\f$ the total Compton scattering cross section as defined
    earlier, \f${\rm S}_{11}(x,\theta)\f$ and \f${\rm S}_{12}(x,\theta)\f$ two of the Compton
    Müller matrix elements, and \f$P_\text{L}\f$ and \f$\gamma\f$ the incoming photon's linear
    polarization degree and angle.

    The marginal distribution for \f$\theta\f$ is determined by the first term of this equation
    (the second term cancels out during the integration over \f$\varphi\f$), which is identical to
    the univariate distribution for \f$\theta\f$ given earlier for the unpolarized case. We can
    thus use the same sampling procedure for \f$\theta\f$.

    Once a random \f$\theta\f$ has been selected, we sample an azimuth angle \f$\varphi\f$ from the
    normalized conditional distribution, again obtained similary as in Peest et al. 2017:

    \f[ \Phi'_\theta(x, \varphi, {\bf{S}}) \propto 1 + \frac{{\rm S}_{12}(x,\theta)} {{\rm
    S}_{11}(x,\theta)} P_\text{L}\cos2(\varphi-\gamma) \f]

    This expression has the same form as the formula derived for Thomson scattering by Peest et al.
    2017, but with different matrix elements \f${\rm S}_{11}(x,\theta)\f$ and \f${\rm
    S}_{12}(x,\theta)\f$. Given a uniform deviate \f$\chi\f$ between 0 and 1, a random
    \f$\varphi\f$ can be obtained from this distribution by solving the equation

    \f[ \chi = \int_0^\varphi \Phi'_\theta(x, \varphi', {\bf{S}})\, d\varphi' =
    \frac{1}{2\pi}\left( \varphi + \frac{{{\rm S}_{12}(\theta, x)}}{{\rm S}_{11}(\theta, x)}
    P_\text{L}\sin\varphi\cos(\varphi-2\gamma)\right), \f]

    which must be inverted for \f$\varphi\f$ numerically. */
class ComptonPhaseFunction
{
    //============= Construction - Setup - Destruction =============

public:
    /** This function caches a pointer to the simulation's random number generator and, if
        polarization is included, pre-calculates some data used for sampling from the relevant
        phase function. The function must be called during setup (i.e. in single-thread mode). If
        \em includePolarization is omitted or set to false, calling the functions implementing the
        polarized phase function will cause undefined behavior. */
    void initialize(Random* random, bool includePolarization = false);

    //======== Absorption =======

    /** This function calculates and returns the Compton scattering cross section for the given
        wavelength, normalized to the (constant) Thomson cross section. */
    double sectionSca(double lambda) const;

    //======== Scattering without polarization =======

private:
    /** This function returns the value of the scattering phase function \f$\Phi(x, \cos\theta)\f$
        for the specified incoming photon energy \f$x\f$ and scattering angle cosine
        \f$\cos\theta\f$, where the phase function is normalized for all \f$x\f$ as \f[\int_{-1}^1
        \Phi(x, \cos\theta) \,\mathrm{d}\cos\theta =2.\f] */
    double phaseFunctionValueForCosine(double x, double costheta) const;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi(x, \cos\theta)\f$ for the specified incoming photon energy \f$x\f$. */
    double generateCosineFromPhaseFunction(double x) const;

    //======== Scattering with polarization =======

private:
    /** This function returns the value of the scattering phase function \f$\Phi(\theta,\phi)\f$
        for the specified incoming photon energy, the specified scattering angles \f$\theta\f$ and
        \f$\phi\f$, and the specified incoming polarization state. The phase function is normalized
        as \f[\int\Phi(\theta,\phi) \,\mathrm{d}\Omega =4\pi.\f]

        For Compton scattering, ... */
    double phaseFunctionValue(double x, double theta, double phi, const StokesVector* sv) const;

    /** This function generates random scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from
        the phase function \f$\Phi(\theta,\phi)\f$ for the specified incoming photon energy and the
        specified incoming polarization state. The results are returned as a pair of numbers in the
        order \f$\theta\f$ and \f$\phi\f$.

        For Compton scattering, ... */
    std::pair<double, double> generateAnglesFromPhaseFunction(double x, const StokesVector* sv) const;

    /** This function applies the Mueller matrix transformation for the specified incoming photon
        energy and the specified scattering angle \f$\theta\f$ to the given polarization state
        (which serves as both input and output for the function).

        For Compton scattering, ... */
    void applyMueller(double x, double theta, StokesVector* sv) const;

    //======== Perform scattering with or without polarization =======

public:
    /** This function calculates the contribution of a Compton scattering event to the peel-off
        photon luminosity and polarization state and determines the adjusted wavelength of the
        outgoing photon packet for the given geometry and incoming wavelength and polarization
        state. The contributions to the Stokes vector components are added to the incoming values
        of the \em I, \em Q, \em U, \em V arguments, and the adjusted wavelength is stored in the
        \em lambda argument. See the description of the MaterialMix::peeloffScattering() function
        for more information. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfk, Direction bfkobs,
                           Direction bfky, const StokesVector* sv) const;

    /** Given the incoming photon packet wavelength, direction  and polarization
        state this function calculates a
        randomly sampled new propagation direction for a Compton scattering event, and determines
        the adjusted wavelength of the outgoing photon packet. The adjusted wavelength is stored in
        the \em lambda argument, and the direction is returned. */
    Direction performScattering(double& lambda, Direction bfk, StokesVector* sv) const;

    //======================== Data Members ========================

private:
    // the simulation's random number generator - initialized by initialize()
    Random* _random{nullptr};
    bool _includePolarization{false};
};

////////////////////////////////////////////////////////////////////

#endif

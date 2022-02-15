/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DIPOLEPHASEFUNCTION
#define DIPOLEPHASEFUNCTION

#include "Array.hpp"
#include "Direction.hpp"
class Random;
class StokesVector;

////////////////////////////////////////////////////////////////////

/** The DipolePhaseFunction helper class represents the dipole scattering phase function, with
    optional support for polarization by scattering. For example, in a wide wavelength range,
    scattering of photons by electrons can be described by elastic and wavelength-independent
    Thomson scattering, which is in turn described using the dipole phase function. The
    corresponding Mueller matrix elements (and hence the phase function) can be expressed
    analytically as a function of just the scattering angle; see Bohren & Huffman (1998) or Wolf
    2003 (Computer Physics Communications, 150, 99–115).

    The implementation of this class is based on the analysis presented by Peest at al. 2017 (A&A,
    601, A92). */
class DipolePhaseFunction
{
    //============= Construction - Setup - Destruction =============

public:
    /** This function caches a pointer to the simulation's random number generator and, if
        polarization is included, pre-calculates some data used for sampling from the corresponding
        phase function. The function must be called during setup (i.e. in single-thread mode). If
        \em includePolarization is omitted or set to false, calling the functions implementing the
        polarized phase function will cause undefined behavior. */
    void initialize(Random* random, bool includePolarization = false);

    //======== Scattering without polarization =======

private:
    /** This function returns the value of the scattering phase function \f$\Phi(\cos\theta)\f$ for
        the specified scattering angle cosine \f$\cos\theta\f$, where the phase function is
        normalized as \f[\int_{-1}^1 \Phi(\cos\theta) \,\mathrm{d}\cos\theta =2.\f]

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero, which yields simply \f[ \Phi(\theta) =
        \frac{3}{4}\,(\cos^2\theta+1). \f] */
    double phaseFunctionValueForCosine(double costheta) const;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi(\cos\theta)\f$.

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero, which yields \f[ \Phi(\theta) \propto \cos^2\theta+1. \f] We can
        thus use the same procedure for sampling \f$\theta\f$ as described for the
        phaseFunctionValue() function. */
    double generateCosineFromPhaseFunction() const;

    //======== Scattering with polarization =======

private:
    /** This function returns the value of the scattering phase function \f$\Phi(\theta,\phi)\f$
        for the specified scattering angles \f$\theta\f$ and \f$\phi\f$, and for the specified
        incoming polarization state. The phase function is normalized as \f[\int\Phi(\theta,\phi)
        \,\mathrm{d}\Omega =4\pi.\f]

        For dipole scattering, we can substitute the analytical expressions for the Mueller
        coefficients into the phase function, which leads to \f[ \Phi(\theta,\phi) =
        \frac{3}{4}\,\left[(\cos^2\theta+1) + P_{\text{L}}\,(\cos^2\theta-1)\cos2(\phi - \gamma)
        \right], \f] where \f$P_\text{L}\f$ is the linear polarization degree and \f$\gamma\f$ the
        polarization angle of the incoming photon. */
    double phaseFunctionValue(double theta, double phi, const StokesVector* sv) const;

    /** This function generates random scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from
        the phase function \f$\Phi(\theta,\phi)\f$ for the specified incoming polarization state.
        The results are returned as a pair of numbers in the order \f$\theta\f$ and \f$\phi\f$.

        For dipole scattering, we sample from the phase function listed for the
        phaseFunctionValue() function using the conditional probability technique. We reduce the
        phase function to the marginal distribution by integrating over \f$\phi\f$, yielding \f[
        \Phi(\theta) \propto \cos^2\theta+1. \f] With the substitution \f$t = \cos\theta\f$, the
        normalized probability distribution to be sampled can be written as \f[ p(t) =
        \frac{3}{8}\,(t^2+1),\text{ with } -1<t<1. \f] The corresponding cumulative distribution is
        \f[ P(t) = \int_{-1}^t \frac{3}{8}\,(t'^2+1)\,\text{d}t' = \frac{1}{8}\,(t^3+3t+4). \f] The
        equation \f${\cal{X}} = P(t)\f$, with \f${\cal{X}}\f$ a uniform deviate, has a single
        real-valued solution, \f[ t = p - \frac{1}{p} \quad\text{with}\quad p = \left( 4{\cal{X}}-2
        +\sqrt{16{\cal{X}}^2 -16{\cal{X}} +5} \right)^{1/3}\f]

        Once we have selected a random scattering angle \f$\theta=\arccos(t)\f$, we sample a random
        azimuthal angle \f$\phi\f$ from the normalized conditional distribution, \f[
        \Phi_\theta(\phi) =\frac{\Phi(\theta,\phi)}{\int_0^{2\pi}
        \Phi(\theta,\phi')\,\text{d}\phi'} =\frac{1}{2\pi}\left(1+
        P_{\text{L}}\,\frac{\cos^2\theta-1}{\cos^2\theta+1}\cos 2(\phi - \gamma)\right). \f]

        This can be done by numerically solving the equation \f[ {\cal{X}}
        =\int_{0}^{\phi}\Phi_{\theta}(\phi')\,\text{d}\phi' =\frac{1}{2\pi} \left( \phi +
        P_{\text{L}}\,\frac{\cos^2\theta-1}{\cos^2\theta+1} \sin\phi \cos(\phi - 2\gamma)\right)
        \f] for \f$\phi\f$, with \f${\cal{X}}\f$ being a new uniform deviate. */
    std::pair<double, double> generateAnglesFromPhaseFunction(const StokesVector* sv) const;

    /** This function applies the Mueller matrix transformation for the specified scattering angle
        \f$\theta\f$ to the given polarization state (which serves as both input and output for the
        function).

        For scattering by spherical particles, the Mueller matrix has only four independent
        coefficients. For dipole scattering, these can be analytically expressed as
        \f[S_{11}=\cos^2\theta+1,\quad S_{12}=\cos^2\theta-1,\quad S_{33}=2\cos\theta,\quad \text{
        and } S_{34}=0.\f] */
    void applyMueller(double theta, StokesVector* sv) const;

    //======== Perform scattering with or without polarization =======

public:
    /** This function calculates the contribution of a dipole scattering event to the peel-off
        photon luminosity and polarization state for the given geometry and polarization state. The
        contributions to the Stokes vector components are added to the incoming values of the \em
        I, \em Q, \em U, \em V arguments. See the description of the
        MaterialMix::peeloffScattering() function for more information. */
    void peeloffScattering(double& I, double& Q, double& U, double& V, Direction bfk, Direction bfkobs, Direction bfky,
                           const StokesVector* sv) const;

    /** Given the incoming photon packet direction and polarization state, this function calculates
        and returns a randomly sampled new propagation direction for a dipole scattering event, and
        if applicable (depending in the polarization flag passed to the constructor), updates the
        polarization state of the photon packet along the way. */
    Direction performScattering(Direction bfk, StokesVector* sv) const;

    //======================== Data Members ========================

private:
    // the simulation's random number generator - initialized during construction
    Random* _random{nullptr};
    bool _includePolarization{false};

    // precalculated discretizations - initialized during construction
    Array _phiv;   // indexed on f
    Array _phi1v;  // indexed on f
    Array _phisv;  // indexed on f
    Array _phicv;  // indexed on f
};

////////////////////////////////////////////////////////////////////

#endif

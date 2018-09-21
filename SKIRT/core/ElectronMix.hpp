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
    Physics Communications, 150, 99–115).

    The implementation of this class is based on the analysis presented by Peest at al. 2017 (A&A,
    601, A92). */
class ElectronMix : public MaterialMix
{
    ITEM_CONCRETE(ElectronMix, MaterialMix, "a population of electrons")

    PROPERTY_BOOL(includePolarization, "include support for polarization")
        ATTRIBUTE_DEFAULT_VALUE(includePolarization, "false")
        ATTRIBUTE_DISPLAYED_IF(includePolarization, "Level2")

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
        \Phi_\lambda(\cos\theta) \,\mathrm{d}\cos\theta =2.\f]

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero, which yields simply \f[ \Phi(\theta) =
        \frac{3}{4}\,(\cos^2\theta+1). \f] */
    double phaseFunctionValueForCosine(double lambda, double costheta) const override;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$.

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero, which yields \f[ \Phi(\theta) \propto \cos^2\theta+1. \f] We can
        thus use the same procedure for sampling \f$\theta\f$ as described for the
        phaseFunctionValue() function. */
    double generateCosineFromPhaseFunction(double lambda) const override;

    //======== Polarization through scattering by spherical particles =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angles \f$\theta\f$ and \f$\phi\f$, and for the specified incoming polarization state. The
        phase function is normalized as \f[\int\Phi_\lambda(\theta,\phi) \,\mathrm{d}\Omega
        =4\pi.\f]

        For Thomson scattering, we can substitute the analytical expressions for the Mueller
        coefficients into the phase function, which leads to \f[ \Phi(\theta,\phi) =
        \frac{3}{4}\,\left[(\cos^2\theta+1) + P_{\text{L}}\,(\cos^2\theta-1)\cos2(\phi - \gamma)
        \right], \f] where \f$P_\text{L}\f$ is the linear polarization degree and \f$\gamma\f$ the
        polarization angle of the incoming photon, and there is no wavelength dependence. */
   double phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const override;

    /** This function generates random scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from
        the phase function \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$, and for the
        specified incoming polarization state. The results are returned as a pair of numbers in the
        order \f$\theta\f$ and \f$\phi\f$.

        For Thomson scattering, we sample from the phase function listed for the
        phaseFunctionValue() function using the conditional probability technique. We reduce the
        phase function to the marginal distribution by integrating over \f$\phi\f$, yielding \f[
        \Phi(\theta) \propto \cos^2\theta+1. \f] With the substitution \f$t = \cos\theta\f$, the
        normalized probability distribution to be sampled can be written as \f[ p(t) =
        \frac{3}{8}\,(t^2+1),\text{ with } -1<t<1. \f] The corresponding cumulative distribution is
        \f[ P(t) = \int_{-1}^t \frac{3}{8}\,(t'^2+1)\,\text{d}t' = \frac{1}{8}\,(t^3+3t+4). \f] The
        equation \f${\cal{X}} = P(t)\f$, with \f${\cal{X}}\f$ a uniform deviate, has a single
        real-valued solution, \f[ t = p - \frac{1}{p} \quad\text{with}\quad p =
        \left( 4{\cal{X}}-2 +\sqrt{16{\cal{X}}^2 -16{\cal{X}} +5} \right)^{1/3}\f]

        Once we have selected a random scattering angle \f$\theta=\arccos(t)\f$, we sample a random
        azimuthal angle \f$\phi\f$ from the normalized conditional distribution, \f[
        \Phi_\theta(\phi) =\frac{\Phi(\theta,\phi)}{\int_0^{2\pi}
        \Phi(\theta,\phi')\,\text{d}\phi'} =\frac{1}{2\pi}\left(1+
        P_{\text{L}}\,\frac{\cos^2\theta-1}{\cos^2\theta+1}\cos 2(\phi - \gamma)\right). \f]

        This can again be done through numerical inversion, by solving the equation \f[ {\cal{X}}
        =\int_{0}^{\phi}\Phi_{\theta}(\phi')\,\text{d}\phi' =\frac{1}{2\pi} \left( \phi +
        P_{\text{L}}\,\frac{\cos^2\theta-1}{\cos^2\theta+1} \sin\phi \cos(\phi - 2\gamma)\right)
        \f] for \f$\phi\f$, with \f${\cal{X}}\f$ being a new uniform deviate. */
    std::pair<double,double> generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const override;

    /** This function applies the Mueller matrix transformation for the specified wavelength
        \f$\lambda\f$ and scattering angle \f$\theta\f$ to the given polarization state (which
        serves as both input and output for the function).

        For scattering by spherical particles, the Mueller matrix has only four independent
        coefficients. For Thomson scattering, these can be analytically expressed as
        \f[S_{11}=\cos^2\theta+1,\quad S_{12}=\cos^2\theta-1,\quad S_{33}=2\cos\theta,\quad \text{
        and } S_{34}=0,\f] i.e. without dependence on wavelength. */
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

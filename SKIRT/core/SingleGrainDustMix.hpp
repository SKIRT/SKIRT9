/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEGRAINDUSTMIX_HPP
#define SINGLEGRAINDUSTMIX_HPP

#include "MaterialMix.hpp"
#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

/** SingleGrainDustMix is an abstract class implementing a dust mix described by a single
    representative grain, with or without support for polarization by scattering. This base class
    includes the implementations of the required functions. Subclasses must merely provide the
    names of the relevant resource files and implement the scatteringMode() function (the
    implementation in this class auto-adapts to the returned value).

    The implementation of this class is based on the analysis presented by Peest at al. 2017 (A&A,
    601, A92). */
class SingleGrainDustMix : public MaterialMix
{
    ITEM_ABSTRACT(SingleGrainDustMix, MaterialMix, "a dust mix described by a single representative grain")
        ATTRIBUTE_TYPE_INSERT(SingleGrainDustMix, "Dust")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function opens the stored table resource(s) tabulating the representative grain
        properties and caches some pre-calculated arrays, depending on the value returned by the
        scatteringMode() function implemented in each subclass. It invokes the resourceNameXxx()
        functions provided by each subclass to obtain the appropriate resource names. */
    void setupSelfBefore() override;

    /** This function must be implemented in a subclass to return the name of the stored table
        resource tabulating the basic optical properties (cross sections and asymmetry parameter)
        as a function of wavelength. The asymmetry parameter values are used only with the
        HenyeyGreenstein scattering mode. */
    virtual string resourceNameForOpticalProps() const = 0;

    /** This function must be implemented in a subclass to return the name of the stored table
        resource tabulating the coefficients of the Mueller matrix as a function of wavelength and
        scattering angle. This function is invoked for the MaterialPhaseFunction scattering mode
        (to obtain \f$S_{11}\f$ ) and for the SphericalPolarization (to obtain all four
        \f$S_{xx}\f$ coefficients). The default implementation in this base class returns the empty
        string, which is acceptable only with the HenyeyGreenstein scattering mode. */
    virtual string resourceNameForMuellerMatrix() const;

    //======== Functionality levels =======

public:
    /** This function returns the fundamental material type represented by this material mix, in
        other words it returns MaterialType::Dust. See the documentation of the MaterialMix class
        for more information. */
    MaterialType materialType() const override;

    //======== Basic material properties =======

public:
    /** This function returns the dust mass per hydrogen atom \f$\mu\f$ for this material. */
    double mass() const override;

    /** This function returns the absorption cross section per hydrogen atom
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sectionAbsSelf(double lambda) const override;

    /** This function returns the scattering cross section per hydrogen atom
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sectionScaSelf(double lambda) const override;

    /** This function returns the scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$, or if this value is unkown, it
        returns zero (corresponding to isotropic scattering). */
    double asymmpar(double lambda) const override;

    //======== Scattering with material phase function =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angle cosine \f$\cos\theta\f$, where the phase function is normalized as \f[\int_{-1}^1
        \Phi_\lambda(\cos\theta) \,\mathrm{d}\cos\theta =2.\f]

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero. The simplified function uses only the first Mueller matrix
        coefficient \f$S_{11}\f$. */
    virtual double phaseFunctionValueForCosine(double lambda, double costheta) const override;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$.

        The phase function for unpolarized radiation is obtained from the general polarized case
        described for the phaseFunctionValue() function by setting the linear polarization degree
        \f$P_\text{L}\f$ to zero. The simplified function no longer depends on \f$\phi\f$, and uses
        only the first Mueller matrix coefficient \f$S_{11}\f$. To sample the scattering angle
        \f$\theta\f$ from this phase function, we follow the same procedure as described for the
        generateAnglesFromPhaseFunction() function, with \f$P_\text{L}=0\f$. */
    virtual double generateCosineFromPhaseFunction(double lambda) const override;

    //======== Polarization through scattering by spherical particles =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angles \f$\theta\f$ and \f$\phi\f$, and for the specified incoming polarization state. The
        phase function is normalized as \f[\int\Phi_\lambda(\theta,\phi) \,\mathrm{d}\Omega
        =4\pi.\f]

        The phase function for scattering by spherical grains can be written as \f[
        \Phi(\theta,\phi) = N\,S_{11} \left( 1 + P_{\text{L}}\,\frac{S_{12}}{S_{11}}\cos2(\phi -
        \gamma) \right) \f] with \f[ N=\frac{1}{2\pi\int_0^\pi S_{11}\sin\theta\, \text{d}\theta},
        \f] where \f$P_\text{L}\f$ is the linear polarization degree and \f$\gamma\f$ the
        polarization angle of the incoming photon, and where the Mueller matrix coefficients
        \f$S_{xx}\f$ depend on both the photon wavelength \f$\lambda\f$ and the scattering angle
        \f$\theta\f$. As a result, the normalization constant cannot be pre-calculated. */
   double phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const override;

    /** This function generates random scattering angles \f$\theta\f$ and \f$\phi\f$ sampled from
        the phase function \f$\Phi_\lambda(\theta,\phi)\f$ at wavelength \f$\lambda\f$, and for the
        specified incoming polarization state. The results are returned as a pair of numbers in the
        order \f$\theta\f$ and \f$\phi\f$.

        For scattering by spherical grains, we sample from the phase function listed for the
        phaseFunctionValue() function using the conditional probability technique. We reduce the
        phase function to the marginal distribution \f$\Phi(\theta)\f$, \f[ \Phi(\theta)
        =\int_0^{2\pi} \Phi(\theta,\phi)\,\text{d}\phi =2\pi\ N\,S_{11} =\frac{S_{11}}{\int_0^\pi
        S_{11}\sin\theta'\, \text{d}\theta'} \ . \f] We sample a random \f$\theta\f$ value from
        this distribution through numerical inversion, that is to say, by solving the equation
        \f[\label{eq:numInvTheta} {\cal{X}} =\frac{\int_0^\theta
        S_{11}\sin\theta'\,\text{d}\theta'}{\int_0^\pi S_{11}\sin\theta'\, \text{d}\theta'} \f] for
        \f$\theta\f$, where \f${\cal{X}}\f$ is a uniform deviate.

        Once we have selected a random scattering angle \f$\theta\f$, we sample a random azimuthal
        angle \f$\phi\f$ from the normalized conditional distribution, \f[ \Phi_\theta(\phi)
        =\frac{\Phi(\theta,\phi)}{\int_0^{2\pi} \Phi(\theta,\phi')\,\text{d}\phi'}
        =\frac{1}{2\pi}\left(1+ P_{\text{L}}\,\frac{S_{12}}{S_{11}}\cos 2(\phi - \gamma)\right),
        \f] where the ratio \f$S_{12}/S_{11}\f$ depends on the scattering angle \f$\theta\f$ in
        addition to wavelength.

        This can again be done through numerical inversion, by solving the equation \f[ {\cal{X}}
        =\int_{0}^{\phi}\Phi_{\theta}(\phi')\,\text{d}\phi' =\frac{1}{2\pi} \left( \phi +
        P_{\text{L}}\,\frac{S_{12}}{S_{11}} \sin\phi \cos(\phi - 2\gamma)\right) \f] for
        \f$\phi\f$, with \f${\cal{X}}\f$ being a new uniform deviate. */
    std::pair<double,double> generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const override;

    /** This function applies the Mueller matrix transformation for the specified wavelength
        \f$\lambda\f$ and scattering angle \f$\theta\f$ to the given polarization state (which
        serves as both input and output for the function).

        For scattering by spherical grains, the Mueller matrix has only four independent
        coefficients, namely \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{33}\f$, and \f$S_{34}\f$, which
        depend on both the photon wavelength \f$\lambda\f$ and the scattering angle \f$\theta\f$.
        These coefficients are obtained from the tables loaded during setup (and specified by name
        in each subclass. */
    void applyMueller(double lambda, double theta, StokesVector* sv) const override;

    //======================== Data Members ========================

private:
    // basic properties - initialized during setup
    StoredTable<1> _sigmaabs;       // indexed on lambda
    StoredTable<1> _sigmasca;       // indexed on lambda
    StoredTable<1> _asymmpar;       // indexed on lambda
    double _mu{0};

    // Mueller matrix - initialized during setup
    StoredTable<2> _S11;            // indexed on lambda, theta
    StoredTable<2> _S12;            // indexed on lambda, theta
    StoredTable<2> _S33;            // indexed on lambda, theta
    StoredTable<2> _S34;            // indexed on lambda, theta

    // precalculated discretizations - initialized during setup
    Array _thetav;                  // indexed on t
    Array _thetasv;                 // indexed on t
    Array _phiv;                    // indexed on f
    Array _phi1v;                   // indexed on f
    Array _phisv;                   // indexed on f
    Array _phicv;                   // indexed on f
};

////////////////////////////////////////////////////////////////////

#endif

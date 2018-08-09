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
    implementation in this class auto-adapts to the returned value). */
class SingleGrainDustMix : public MaterialMix
{
    ITEM_ABSTRACT(SingleGrainDustMix, MaterialMix, "a dust mix described by a single representative grain")
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
        resource tabulating the elements of the Mueller matrix as a function of wavelength and
        scattering angle. This function is invoked for the MaterialPhaseFunction scattering mode
        (to obtain S11) and for the SphericalPolarization (to obtain all four Sxx coefficients).
        The default implementation in this base class returns the empty string, which is acceptable
        only with the HenyeyGreenstein scattering mode. */
    virtual string resourceNameForMuellerMatrix() const;

    //======== Functionality levels =======

public:
    /** This function returns the fundamental material type represented by this material mix, in
        other words it returns MaterialType::Dust. See the documentation of the MaterialMix class
        for more information. */
    MaterialType materialType() const override;

    //======== Basic material properties =======

public:
    /** This function returns the absorption cross section per hydrogen atom
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per hydrogen atom
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sectionSca(double lambda) const override;

    /** This function returns the total extinction cross section per hydrogen atom
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ of the dust mix at
        wavelength \f$\lambda\f$. */
    double albedo(double lambda) const override;

    /** This function returns the scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$, or if this value is unkown, it
        returns zero (corresponding to isotropic scattering). */
    double asymmpar(double lambda) const override;

    //======== Scattering with material phase function =======

public:
    /** This function returns the value of the scattering phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$ for the specified scattering
        angle cosine \f$\cos\theta\f$, where the phase function is normalized as \f[\int_{-1}^1
        \Phi_\lambda(\cos\theta) \,\mathrm{d}\cos\theta =2.\f] */
    virtual double phaseFunctionValueForCosine(double lambda, double costheta) const override;

    /** This function generates a random scattering angle cosine sampled from the phase function
        \f$\Phi_\lambda(\cos\theta)\f$ at wavelength \f$\lambda\f$. */
    virtual double generateCosineFromPhaseFunction(double lambda) const override;

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

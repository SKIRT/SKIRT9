/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIGRAINDUSTMIX_HPP
#define MULTIGRAINDUSTMIX_HPP

#include "DustMix.hpp"
#include "GrainPopulation.hpp"
class GrainComposition;
class GrainSizeDistribution;

////////////////////////////////////////////////////////////////////

/** MultiGrainDustMix is an abstract class implementing a dust mix described by one or more grain
    populations, each with their own grain composition and size distribution, and with or without
    support for polarization by scattering. The class offers facilities to its subclasses to add
    dust grain populations to the dust mix during initial setup. Subsequently, the class uses this
    list of grain populations to calculate the optical properties requested by the DustMix base
    class, and to implement the additional functionality offered by multi-grain-population dust
    mixes, such as providing the information needed for calculating emission from stochastically
    heated dust grains.

    <b>Calculating representative grain properties</b>

    The getOpticalProperties() function calculates the basic representative grain properties
    expected by the DustMix base class from the following information obtained from the dust grain
    populations added by the subclass: the absorption efficiencies \f$Q^{\text{abs}}(\lambda,a)\f$,
    the scattering efficiencies \f$Q^{\text{sca}}(\lambda,a)\f$, the scattering phase function
    asymmetry parameter \f$g(\lambda,a)\f$, the Mueller matrix coefficients
    \f$S^\text{xx}(\lambda,a,\theta)\f$, the bulk density \f$\rho_{\text{bulk}}\f$ of the grain
    material, and the properly normalized grain size distribution per hydrogen atom
    \f$\Omega(a)=(\frac{\text{d}n_\text{D}}{\text{d}a})/n_\text{H}\f$ in the range
    \f$[a_\text{min},a_\text{max}]\f$.

    The properties are calculated for each wavelength \f$\lambda_\ell\f$ in the wavelength grid
    specified by the DustMix class. The values are obtained by integrating over the grain size
    distribution \f$\Omega(a)\f$ using a builtin logarithmic grid and accumulating over all grain
    populations \f$c\f$, using the formulas listed below.

    The absorption and scattering cross sections per hydrogen atom
    \f$\varsigma_{\ell}^{\text{abs}}\f$ and \f$\varsigma_{\ell}^{\text{abs}}\f$ for the
    \f$\ell\f$'th wavelength are calculated using \f[ \varsigma_{\ell}^{\text{abs}} = \sum_c
    \int_{a_{\text{min},c}}^{a_{\text{max},c}} \Omega_c(a)\, Q^{\text{abs}}_c(\lambda_\ell,a)\, \pi
    a^2\, {\text{d}}a \f] and \f[ \varsigma_{\ell}^{\text{sca}} = \sum_c
    \int_{a_{\text{min},c}}^{a_{\text{max},c}} \Omega_c(a)\, Q^{\text{sca}}_c(\lambda_\ell,a)\, \pi
    a^2\, {\text{d}}a. \f]

    The Mueller matrix coefficients provided by the grain population are assumed to be expressed as
    a cross section (in arbitrary units). They are thus integrated over the size distribution
    without again multiplying by the grain cross section, i.e. using \f[
    S^\text{xx}_{\ell,\text{t}} = \sum_c \int_{a_{\text{min},c}}^{a_{\text{max},c}} \Omega_c(a)\,
    S^\text{xx}_c(\lambda_\ell,a,\theta_\text{t})\, {\text{d}}a \f]

    The representative asymmetry parameter \f$g_{\ell}\f$ is averaged over the scattering cross
    section and thus calculated using \f[ g_{\ell} = \frac{1}{\varsigma_{\ell}^{\text{sca}}} \sum_c
    \int_{a_{\text{min},c}}^{a_{\text{max},c}} \Omega_c(a)\, g_c(\lambda_\ell,a)\,
    Q^{\text{sca}}_c(\lambda_\ell,a)\, \pi a^2\, {\text{d}}a. \f]

    The dust mass per hydrogen atom \f$\mu\f$ is calculated by integrating the bulk density over
    the size distribution, \f[ \mu = \sum_c \int_{a_{\text{min},c}}^{a_{\text{max},c}} \Omega_c(a)\,
    \rho_{\text{bulk},c}\, \frac{4\pi}{3}\, a^3\, {\text{d}}a. \f]

    <b>Exposing multiple grain populations</b>

    TO DO: add further documentation.

    <b>Supporting stochastic heating</b>

    TO DO: add further documentation. */
class MultiGrainDustMix : public DustMix
{
    ITEM_ABSTRACT(MultiGrainDustMix, DustMix, "a dust mix with one or more grain populations")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

    //------------- To be invoked by subclasses ------------

protected:
    /** This function adds the specified grain population to the dust mix. The receiving dust mix
        object retains a pointer to the specified GrainPopulation instance for later reference, but
        does not take ownership. The caller must ensure that the GrainPopulation instance lives at
        least as long as the dust mix, and that it is eventually destroyed (at the same time as or
        later than the dust mix). */
    void addPopulation(const GrainPopulation* population);

    /** This function adds the a grain population to the dust mix specified by its constituent
        components. The function creates a new GrainPopulation instance, passing its own function
        arguments to the GrainPopulation constructor. The receiving dust mix object claims
        ownership of the new GrainPopulation instance. The caller must guarantee that the lifetime
        of the specified composition and size distribution objects is as least as long as the
        lifetime of the receiving dust mix.

        Refer to the description of the GrainPopulation constructor for more information on the
        arguments of this function. */
    void addPopulation(GrainComposition* composition, GrainSizeDistribution* sizeDistribution,
                       int numSizes, GrainPopulation::NormalizationType normType, double normValue);

    //------------- Invoked by the DustMix base class ------------

protected:
    /** This function is invoked by the DustMix base class to obtain the representative grain
        optical properties for this dust mix. The first two arguments respectively specify the
        wavelength grid and (if applicable) the scattering angle grid on which the properties must
        be tabulated. The output arrays and tables will already have the appropriate size
        (corresponding to the input wavelength grids) when the function gets called.

        For this class, this function integrates the optical properties over the grain size
        distribution for each of the grain populations added by a subclass as described in the
        class header, and stores the results into the corresponding output arrays. Also, the
        function returns the dust mass per hydrogen atom for the dust mix.

        For the HenyeyGreenstein scattering mode, the Mueller matric tables remain untouched. For
        the MaterialPhaseFunction scattering mode, the function fills only the first table and
        leaves the other tables untouched. For the SphericalPolarization scattering mode, the
        function fills all four tables. */
    double getOpticalProperties(const Array& lambdav, const Array& thetav,
                                Array& sigmaabsv, Array& sigmascav, Array& asymmparv,
                                Table<2>& S11vv, Table<2>& S12vv, Table<2>& S33vv, Table<2>& S34vv) const override;

    /** This function returns the scattering mode supported by this material mix as determined by
        the grain populations added by a subclass. If none of the populations offer Mueller matrix
        coefficients, the function returns the HenyeyGreenstein scattering mode. If all populations
        offer Mueller matrix coefficients, the function returns the SphericalPolarization
        scattering mode. Because there is no configuration option to choose between
        MaterialPhaseFunction or SphericalPolarization in this case, the current implementation
        never returns the MaterialPhaseFunction scattering mode.

        All grain populations should have the same level of Mueller matrix support. If this is not
        the case, this function throws a fatal error. */
    ScatteringMode scatteringMode() const override;

    //======================== Data Members ========================

private:
    // list created by addPopulation()
    vector<const GrainPopulation*> _populations;
};

////////////////////////////////////////////////////////////////////

#endif

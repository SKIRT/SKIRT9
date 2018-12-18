/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MULTIGRAINDUSTMIX_HPP
#define MULTIGRAINDUSTMIX_HPP

#include "DustMix.hpp"
#include "GrainPopulation.hpp"
#include "Range.hpp"
#include "StoredTable.hpp"
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

    All grain populations must have a level of Mueller matrix support that matches the value
    returned by the scatteringMode() function (which may be overridden in a subclass).
    Specifically, for the HenyeyGreenstein scattering mode, \em none of the grain populations
    should offer a Mueller matrix. For the MaterialPhaseFunction and SphericalPolarization
    scattering modes, \em all of the grain populations should offer a Mueller matrix. If this is
    not the case, a fatal error will result.

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
    the size distribution, \f[ \mu = \sum_c \int_{a_{\text{min},c}}^{a_{\text{max},c}}
    \Omega_c(a)\, \rho_{\text{bulk},c}\, \frac{4\pi}{3}\, a^3\, {\text{d}}a. \f]

    <b>Supporting dust emission calculations</b>

    The representative grain properties described above and offered by the public MaterialMix
    interface supported by this class are insufficient to accurately calculate dust emission
    spectra for the dust mixture. This is so because the emission spectrum is a nonlinear function
    of (among many other things) the grain size, and thus a single grain cannot accurately
    represent a population with a (potentialy large) range of grain sizes. Furthermore, smaller
    dust grains are often not in local thermal equilibrium, and instead are heated stochastically
    by individual photon absorption events. Modeling emission for these grains involves a
    temperature probability distribution rather than just an equilibrium temperature. The
    calculation needs calorimetric properties of the grain material in addition to optical
    properties.

    It is numerically intractable to handle every possible grain size seperately. Instead, the
    MultiGrainDustMix class discretizes the grain size distribution for each type of grain material
    into a number of consecutive size bins (on a logarithmic scale), and calculates the optical and
    calorimetric properties of a representative grain for each of these bins. The number of bins
    for each type of grain material can be configured by the user. A larger number of bins improves
    the accuracy of the dust emission spectra. On the other hand, the calculation time scales
    roughly linearly with the number of bins.

    The MultiGrainDustMix class offers a set of public functions on top of the regular MaterialMix
    public interface to expose information about the individual grain populations \f$c\f$ and the
    size-discretized bins \f$b\f$ as described above. For now, these functions are declared and
    implemented here. Over time they might be declared as part of a seperate interface, so that
    they might also be implemented by other MaterialMix subclasses (such as perhaps a material mix
    including both hydrogen gas and dust).

    As a result, the MultiGrainDustMix class supports the folowing methods for calculating dust
    emissivity with varying levels of accuracy:

    - using a single representative grain and assuming local thermal equilibrium

    - using multiple representative grains, one for each composition/size bin, still assuming local
    thermal equilibrium for each bin

    - using multiple representative grains without assuming local thermal equilibrium, i.e.
    calculating a temperature probability distribution to take into account stochastically heated
    grains.

 */
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
                                Table<2>& S11vv, Table<2>& S12vv, Table<2>& S33vv, Table<2>& S34vv) override;

    /** This function is invoked by the DustMix base class to precalculate additional dust
        properties that are offered through functions outside of the DustMix interface. The
        argument specifies the wavelength grid on which the properties must be tabulated (i.e. the
        same grid as passed to the getOpticalProperties() function. The function returns the
        number of memory bytes allocated to store extra propertries.

        The function discretizes the grain size distribution for each grain population added to
        this dust mix into a number of consecutive size bins (on a logarithmic scale), and
        calculates the relevant optical and calorimetric properties of a representative grain for
        each of these bins. The number of bins for each type of grain material can be configured by
        the user.

        The bins for the various populations in the dust mix are placed in a single sequence, in
        the same overall order as the populations addressed by the \f$c\f$ indices. Within each
        population, the bins are listed in order of increasing grain size. */
    size_t initializeExtraProperties(const Array& lambdav) override;

    //------------- Cleanup ------------

public:
    /** The destructor destructs any enthalpy stored tables that were opened during setup. */
    ~MultiGrainDustMix();

    //=============== Exposing multiple grain populations ==============

public:
    /** This function returns the number of dust grain populations (with indices \f$c\f$) added to
        this dust mix by the subclass. Each grain population represents the combination of a grain
        composition, providing the optical and calorimetric properties of the grain material, and a
        grain size distribution with some normalization to specify the the amount of dust contained
        in the population. No grain size discretization has been applied to these populations. */
    int numPopulations() const;

    /** This function returns a brief human-readable identifier for the type of grain material
        represented by the population with index \f$c\f$. The identifier does not contain white
        space. */
    string populationGrainType(int c) const;

    /** This function returns the minimum and maximum grain sizes \f$a_{\text{min},c},
        a_{\text{max},c}\f$ for the population with index \f$c\f$. */
    Range populationSizeRange(int c) const;

    /** This function returns the dust mass \f$\mu_c\f$ per hydrogen atom for the population with
        index \f$c\f$. */
    double populationMass(int c) const;

    /** This function returns the number of dust grain size bins (with indices \f$b\f$) in this
        dust mix. During setup, the MultiGrainDustMix class discretizes the grain size distribution
        for each grain population added to this dust mix into a number of consecutive size bins (on
        a logarithmic scale), and calculates the relevant optical and calorimetric properties of a
        representative grain for each of these bins. The number of bins for each type of grain
        material can be configured by the user.

        The bins for the various populations in the dust mix are placed in a single sequence, in
        the same overall order as the populations addressed by the \f$c\f$ indices. Within each
        population, the bins are listed in order of increasing grain size. For example, a dust mix
        with a silicate population discretized into 7 grain size bins, a graphite population with 6
        size bins, and a PAH population with 5 bins, will have a total of 18 bins, consecutively
        addressed by the \f$b\f$ indices. */
    int numBins() const;

    /** This function returns the equilibrium temperature \f$T_{\text{eq},b}\f$ (assuming LTE
        conditions) for the representative grain of the bin with index \f$b\f$ when it would be
        embedded in the radiation field specified by the mean intensities \f$(J_\lambda)_\ell\f$,
        which must be discretized on the simulation's radiation field wavelength grid as returned
        by the Configuration::radiationFieldWLG() function. */
    double binEquilibriumTemperature(int b, const Array& Jv) const;

    /** This function returns the mean mass of a dust grain for the bin with index \f$b\f$. */
    double binMeanMass(int b) const;

    /** This function returns the absorption cross section per hydrogen atom
        \f$\varsigma^{\text{abs}}_{\lambda,b}\f$ at wavelength \f$\lambda\f$ for the representative
        grain of the bin with index \f$b\f$. */
    double binSectionAbs(int b, double lambda) const;

    /** This function returns the enthalpy at temperature \f$T\f$ for the representative grain of
        the bin with index \f$b\f$. The enthalpy is equivalent to the internal energy of the dust
        grain, using an arbitrary zero point. It is obtained by multiplying the specific enthalpy
        of the appropriate grain composition (at the specified temperature) by the mean mass of a
        dust grain in the requested bin. If the specified temperature lies outside of the
        internally defined grid, the enthalpy value at the nearest border is used instead. */
    double binEnthalpy(int b, double T) const;

    /** This function returns the largest temperature for which this dust mix can provide
        meaningful enthalpy data, for any of the dust grain populations in the mix. */
    double maxEnthalpyTemperature() const;

    //======================== Data Members ========================

private:
    // list created by addPopulation()
    vector<const GrainPopulation*> _populations;

    // data members with information per population -- initialized by getOpticalProperties()
    vector<double> _mupopv;             // mass per hydrogen atom for population - indexed on c
    vector<double> _normv;              // size distribution normalization for population - indexed on c

    // data members with info per size bin for multigrain emissivity -- initialized by initializeExtraProperties()
    EquilibriumDustTemperatureCalculator _tempCalcv; // equilibrium temperature info - indexed on b

    // data members with info per size bin for stochastic emissivity -- initialized by initializeExtraProperties()
    vector<int> _btocv;     // mapping from b index (bins) to c index (corresponding population)
    Array _massv;           // mean mass of a grain - indexed on b
    Table<2> _sigmaabsvv;   // absorption cross sections - indexed on b,ell
    vector<StoredTable<1>*> _enthalpyv; // enthalpy stored table for population - indexed on c
};

////////////////////////////////////////////////////////////////////

#endif

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
    support for polarization by scattering.

    TO DO: add further documentation.
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
        distribution for each of the grain populations added by a subclass, and stores the results
        into the corresponding output arrays. Also, the function returns the dust mass per hydrogen
        atom for the dust mix.

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

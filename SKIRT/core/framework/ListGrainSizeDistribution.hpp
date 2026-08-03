/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTGRAINSIZEDISTRIBUTION_HPP
#define LISTGRAINSIZEDISTRIBUTION_HPP

#include "Array.hpp"
#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** ListGrainSizeDistribution represents a grain size distribution that is fully specified inside
    the configuration file (i.e. without referring to an input file). It is intended for use in
    cases where there are just a few data points in the size distribution, but nothing keeps the
    user from specifying a long list. The grain sizes \f$a\f$ must listed be in increasing order,
    and the size distribution values \f$\text{dnda} \propto \frac{\text{d}n_\text{D}}{\text{d}a}\f$
    must be specified corresponding to each grain size. Where needed, the tabulated values are
    interpolated logarithmically on both axes. Outside of the specified grain size range, the
    number of dust grains is considered to be zero.

    Any required normalization can and should be applied externally to this class in the
    configuration of the grain population (see GrainSizeDistribution and GrainPopulation). */
class ListGrainSizeDistribution : public GrainSizeDistribution
{
    ITEM_CONCRETE(ListGrainSizeDistribution, GrainSizeDistribution,
                  "a dust grain size distribution specified inside the configuration file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ListGrainSizeDistribution, "Level2")

        PROPERTY_DOUBLE_LIST(sizes, "the grain sizes at which to specify the size distribution")
        ATTRIBUTE_QUANTITY(sizes, "grainsize")
        ATTRIBUTE_MIN_VALUE(sizes, "[1 Angstrom")
        ATTRIBUTE_MAX_VALUE(sizes, "1 mm]")

        PROPERTY_DOUBLE_LIST(sizeDistributionValues, "the size distribution values at each of the given grain sizes")
        ATTRIBUTE_QUANTITY(sizeDistributionValues, "pergrainsize")
        ATTRIBUTE_MIN_VALUE(sizeDistributionValues, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the number of configured values. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the minimum grain size \f$a_\text{min}\f$, i.e. the first grain size
        in the configured list. */
    double amin() const override;

    /** This function returns the maximum grain size \f$a_\text{max}\f$, i.e. the last grain size
        in the configured list. */
    double amax() const override;

    /** This function returns the value of the distribution \f$\text{dnda} \propto
        \frac{\text{d}n_\text{D}}{\text{d}a}\f$ for a given grain size \f$a\f$, log-log
        interpolated from the configured table. If \f$a<a_\text{min}\f$ or \f$a>a_\text{max}\f$ the
        function returns zero. */
    double dnda(double a) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _av;     // grain sizes
    Array _dndav;  // size distribution values
};

////////////////////////////////////////////////////////////////////

#endif

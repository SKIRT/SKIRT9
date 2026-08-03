/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RANGEGRAINSIZEDISTRIBUTION_HPP
#define RANGEGRAINSIZEDISTRIBUTION_HPP

#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** RangeGrainSizeDistribution is an abstract class that represents a grain size distribution with
    a configurable grain size range. Specifically, it manages the attributes \f$a_\text{min}\f$ and
    \f$a_\text{max}\f$, which determine the range of the distribution. The
    RangeGrainSizeDistribution class consequently implements the functions amin() and amax(), while
    it still expects each subclass to provide the actual distribution function by implementing the
    dnda() function. */
class RangeGrainSizeDistribution : public GrainSizeDistribution
{
    ITEM_ABSTRACT(RangeGrainSizeDistribution, GrainSizeDistribution,
                  "a dust grain size distribution with a configurable size range")

        PROPERTY_DOUBLE(minSize, "the minimum grain size for this distribution")
        ATTRIBUTE_QUANTITY(minSize, "grainsize")
        ATTRIBUTE_MIN_VALUE(minSize, "[1 Angstrom")
        ATTRIBUTE_MAX_VALUE(minSize, "1 mm]")

        PROPERTY_DOUBLE(maxSize, "the maximum grain size for this distribution")
        ATTRIBUTE_QUANTITY(maxSize, "grainsize")
        ATTRIBUTE_MIN_VALUE(maxSize, "[1 Angstrom")
        ATTRIBUTE_MAX_VALUE(maxSize, "1 mm]")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This constructor simply initializes \f$a_\text{min}\f$ and \f$a_\text{max}\f$ with the
        specified values. It can be invoked from subclass constructors to set the values of these
        properties programmatically. */
    explicit RangeGrainSizeDistribution(double minSize, double maxSize);

    /** This function verifies the property values. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the configured minimum grain size \f$a_\text{min}\f$. */
    double amin() const override;

    /** This function returns the configured maximum grain size \f$a_\text{max}\f$. */
    double amax() const override;
};

////////////////////////////////////////////////////////////////////

#endif

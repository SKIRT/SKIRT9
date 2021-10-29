/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SAMPLINGOPTIONS_HPP
#define SAMPLINGOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The SamplingOptions class simply offers a number of configuration options related to sampling
    the media properties for the spatial grid.

    A separate number of spatial samples per grid cell can be configured for the density and for
    all of the the other medium properties (as a group). In each case, if the specified number of
    samples is one, the property is sampled just at the central position of the cell. If the
    specified number is larger than one, the property is sampled at the given number of positions,
    randomly selected from a uniform distribution within the cell volume, and the average value of
    these random samples is used. The average is density-weigthed except for the magnetic field
    (and density itself).

    The medium system maintains at most a single bulk velocity value per spatial cell. If multiple
    medium components specify a bulk velocity, the values sampled from these components are
    aggregated using one of three possible policies, as requested by the user through the \em
    aggregateVelocity property:

    - Average: use the density-weighted average; missing values are taken to be zero.

    - Maximum: use the velocity vector with largest magnitude; missing values are taken to be zero.

    - First: use the vector of the first medium component for which one is available.

    Similarly, the medium system maintains at most a single magnetic field vector per spatial cell.
    However, because the configuration can contain at most one medium component that specifies a
    magnetic field, there is no need for aggregation over multiple components. */
class SamplingOptions : public SimulationItem
{
    /** The enumeration type defining a policy for aggregating (in each spatial cell) a single bulk
        velocity from the values in multiple medium components. */
    ENUM_DEF(AggregatePolicy, Average, Maximum, First)
        ENUM_VAL(AggregatePolicy, Average, "Use the density-weighted average; missing values are taken to be zero")
        ENUM_VAL(AggregatePolicy, Maximum, "Use the vector with largest magnitude; missing values are taken to be zero")
        ENUM_VAL(AggregatePolicy, First, "Use the vector of the first medium component for which one is available")
    ENUM_END()

    ITEM_CONCRETE(SamplingOptions, SimulationItem, "a set of options related to media sampling for the spatial grid")

        PROPERTY_INT(numDensitySamples, "the number of random density samples for determining spatial cell mass")
        ATTRIBUTE_MIN_VALUE(numDensitySamples, "1")
        ATTRIBUTE_MAX_VALUE(numDensitySamples, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numDensitySamples, "100")
        ATTRIBUTE_DISPLAYED_IF(numDensitySamples, "Level2")

        PROPERTY_INT(numPropertySamples, "the number of random samples for determining other medium properties")
        ATTRIBUTE_MIN_VALUE(numPropertySamples, "1")
        ATTRIBUTE_MAX_VALUE(numPropertySamples, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numPropertySamples, "1")
        ATTRIBUTE_DISPLAYED_IF(numPropertySamples, "Level3")

        PROPERTY_ENUM(aggregateVelocity, AggregatePolicy,
                      "aggregating the bulk velocity from multiple medium components")
        ATTRIBUTE_DEFAULT_VALUE(aggregateVelocity, "Average")
        ATTRIBUTE_RELEVANT_IF(aggregateVelocity, "MediumVelocity")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

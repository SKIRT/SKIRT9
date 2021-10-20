/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SAMPLINGOPTIONS_HPP
#define SAMPLINGOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The SamplingOptions class simply offers a number of configuration options related to sampling
    the media properties for the spatial grid. */
class SamplingOptions : public SimulationItem
{
    /** The enumeration type defining a policy for aggregating (in each spatial cell) a single bulk
        velocity or magnetic field vector from the values in multiple medium components. */
    ENUM_DEF(AggregatePolicy, Average, Maximum, Single)
        ENUM_VAL(AggregatePolicy, Average, "Use the average vector; missing values are taken to be zero")
        ENUM_VAL(AggregatePolicy, Maximum, "Use the vector with largest magnitude; missing values are taken to be zero")
        ENUM_VAL(AggregatePolicy, Single, "Require exactly one medium component to provide the value")
    ENUM_END()

    ITEM_CONCRETE(SamplingOptions, SimulationItem, "a set of options related to media sampling for the spatial grid")

        PROPERTY_INT(numDensitySamples, "the number of random density samples for determining spatial cell mass")
        ATTRIBUTE_MIN_VALUE(numDensitySamples, "10")
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
        ATTRIBUTE_DEFAULT_VALUE(velocityPolicy, "Average")
        ATTRIBUTE_RELEVANT_IF(velocityPolicy, "MediumVelocity")

        PROPERTY_ENUM(aggregateMagneticField, AggregatePolicy,
                      "aggregating the magnetic field from multiple medium components")
        ATTRIBUTE_DEFAULT_VALUE(magneticFieldPolicy, "Single")
        ATTRIBUTE_RELEVANT_IF(magneticFieldPolicy, "MagneticField")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

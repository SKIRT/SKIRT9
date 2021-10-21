/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITERATIONOPTIONS_HPP
#define ITERATIONOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The IterationOptions class simply offers a number of options for configuring iterations during
    primary and/or secondary emission. These options are relevant only when the simulation has a
    dynamic medium state and/or a dynamic secondary emission. */
class IterationOptions : public SimulationItem
{
    ITEM_CONCRETE(IterationOptions, SimulationItem,
                  "a set of options for configuring iterations during primary and/or secondary emission")

        PROPERTY_INT(minPrimaryIterations, "the minimum number of iterations during primary emission")
        ATTRIBUTE_MIN_VALUE(minPrimaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(minPrimaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minPrimaryIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minPrimaryIterations, "DynamicState")

        PROPERTY_INT(maxPrimaryIterations, "the maximum number of iterations during primary emission")
        ATTRIBUTE_MIN_VALUE(maxPrimaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxPrimaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxPrimaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxPrimaryIterations, "DynamicState")

        PROPERTY_INT(minSecondaryIterations, "the minimum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(minSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minSecondaryIterations, "DynamicEmission")

        PROPERTY_INT(maxSecondaryIterations, "the maximum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(maxSecondaryIterations, "0")
        ATTRIBUTE_MAX_VALUE(maxSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxSecondaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxSecondaryIterations, "DynamicEmission")

        PROPERTY_BOOL(includePrimaryEmission, "include primary emission in the secondary emission iterations")
        ATTRIBUTE_DEFAULT_VALUE(includePrimaryEmission, "false")
        ATTRIBUTE_RELEVANT_IF(includePrimaryEmission, "DynamicState&DynamicEmission")
        ATTRIBUTE_DISPLAYED_IF(includePrimaryEmission, "Level3")

        PROPERTY_DOUBLE(iterationPacketsMultiplier,
                        "the multiplier on the number of photon packets launched for each iteration")
        ATTRIBUTE_MIN_VALUE(iterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(iterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(iterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(iterationPacketsMultiplier, "DynamicState|DynamicEmission")
        ATTRIBUTE_DISPLAYED_IF(iterationPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

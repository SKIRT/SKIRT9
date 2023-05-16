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
        ATTRIBUTE_RELEVANT_IF(minPrimaryIterations, "IteratePrimary")

        PROPERTY_INT(maxPrimaryIterations, "the maximum number of iterations during primary emission")
        ATTRIBUTE_MIN_VALUE(maxPrimaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxPrimaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxPrimaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxPrimaryIterations, "IteratePrimary")

        PROPERTY_INT(minSecondaryIterations, "the minimum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(minSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minSecondaryIterations, "IterateSecondary")

        PROPERTY_INT(maxSecondaryIterations, "the maximum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(maxSecondaryIterations, "0")
        ATTRIBUTE_MAX_VALUE(maxSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxSecondaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxSecondaryIterations, "IterateSecondary")

        PROPERTY_BOOL(includePrimaryEmission, "include primary emission in the secondary emission iterations")
        ATTRIBUTE_DEFAULT_VALUE(includePrimaryEmission, "false")
        ATTRIBUTE_RELEVANT_IF(includePrimaryEmission, "IterateSecondary")
        ATTRIBUTE_DISPLAYED_IF(includePrimaryEmission, "Level3")

        PROPERTY_DOUBLE(primaryIterationPacketsMultiplier,
                        "the multiplier on the number of photon packets launched for each primary emission iteration")
        ATTRIBUTE_MIN_VALUE(primaryIterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(primaryIterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(primaryIterationPacketsMultiplier, "IteratePrimary|includePrimaryEmission")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationPacketsMultiplier, "Level3")

        PROPERTY_DOUBLE(secondaryIterationPacketsMultiplier,
                        "the multiplier on the number of photon packets launched for each secondary emission iteration")
        ATTRIBUTE_MIN_VALUE(secondaryIterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryIterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryIterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(secondaryIterationPacketsMultiplier, "IterateSecondary")
        ATTRIBUTE_DISPLAYED_IF(secondaryIterationPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

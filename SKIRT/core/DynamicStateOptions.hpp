/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DYNAMICSTATEOPTIONS_HPP
#define DYNAMICSTATEOPTIONS_HPP

#include "DynamicStateRecipe.hpp"
#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The DynamicStateOptions class simply offers a number of configuration options related to the
    capability of dynamically adjusting the medium state by iterating over the radiation field.
    These options are relevant only for simulations that track a panchromatic radiation field with
    a sufficiently broad wavelength range. This includes both extinction-only simulations with the
    flag for storing the radiation field turned on and simulations with secondary emission. On the
    other hand, this excludes oligochromatic simulations and Lya-specific simulations.

    The current implementation supports dynamic medium state iterations during primary emission. If
    required, it probably is fairly straightforward to also support dynamic medium state iterations
    during secondary emission. */
class DynamicStateOptions : public SimulationItem
{
    ITEM_CONCRETE(DynamicStateOptions, SimulationItem,
                  "a set of options for dynamically adjusting the medium state by iterating over the radiation field")

        PROPERTY_ITEM_LIST(recipes, DynamicStateRecipe, "the dynamic medium state recipes")

        PROPERTY_INT(minIterations, "the minimum number of dynamic medium state iterations")
        ATTRIBUTE_MIN_VALUE(minIterations, "1")
        ATTRIBUTE_MAX_VALUE(minIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minIterations, "1")
        ATTRIBUTE_DISPLAYED_IF(minIterations, "Level3")

        PROPERTY_INT(maxIterations, "the maximum number of dynamic medium state iterations")
        ATTRIBUTE_MIN_VALUE(maxIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxIterations, "10")
        ATTRIBUTE_DISPLAYED_IF(maxIterations, "Level3")

        PROPERTY_DOUBLE(
            iterationPacketsMultiplier,
            "the multiplier on the number of photon packets launched for each dynamic medium state iteration")
        ATTRIBUTE_MIN_VALUE(iterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(iterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(iterationPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(iterationPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

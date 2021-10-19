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
    These options are relevant only when this capability has been turned on by the user-configured
    \em iterateMediumState flag in the MonteCarloSimulation class. */
class DynamicStateOptions : public SimulationItem
{
    ITEM_CONCRETE(DynamicStateOptions, SimulationItem, "a set of options for dynamically adjusting the medium state")

        PROPERTY_ITEM_LIST(recipes, DynamicStateRecipe, "the dynamic medium state recipes")
        ATTRIBUTE_RELEVANT_IF(recipes, "DynamicState")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

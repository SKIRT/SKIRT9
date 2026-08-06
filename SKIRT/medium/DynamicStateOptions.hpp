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
    capability of dynamically adjusting the medium state while iterating over primary emission. */
class DynamicStateOptions : public SimulationItem
{
    ITEM_CONCRETE(DynamicStateOptions, SimulationItem, "a set of options for dynamically adjusting the medium state")

        PROPERTY_ITEM_LIST(recipes, DynamicStateRecipe, "the dynamic medium state recipes")
        ATTRIBUTE_RELEVANT_IF(recipes, "IteratePrimary|IterateSecondary")
        ATTRIBUTE_REQUIRED_IF(recipes, "false")
        ATTRIBUTE_INSERT(recipes, "DynamicState")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

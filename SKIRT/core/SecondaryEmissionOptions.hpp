/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SECONDARYEMISSIONOPTIONS_HPP
#define SECONDARYEMISSIONOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The SecondaryEmissionOptions class simply offers a number of configuration options related to
    secondary emission, regardless of the emitting media type. */
class SecondaryEmissionOptions : public SimulationItem
{
    ITEM_CONCRETE(SecondaryEmissionOptions, SimulationItem, "a set of options related to secondary emission")

        PROPERTY_BOOL(storeEmissionRadiationField,
                      "store the radiation field during emission so that it can be probed for output")
        ATTRIBUTE_DEFAULT_VALUE(storeEmissionRadiationField, "false")
        ATTRIBUTE_DISPLAYED_IF(storeEmissionRadiationField, "Level3")

        PROPERTY_DOUBLE(secondaryPacketsMultiplier,
                        "the multiplier on the number of photon packets launched for secondary emission")
        ATTRIBUTE_MIN_VALUE(secondaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(secondaryPacketsMultiplier, "Level3")

        PROPERTY_DOUBLE(spatialBias,
                        "the fraction of secondary photon packets distributed uniformly across spatial cells")
        ATTRIBUTE_MIN_VALUE(spatialBias, "[0")
        ATTRIBUTE_MAX_VALUE(spatialBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(spatialBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(spatialBias, "Level3")

        PROPERTY_DOUBLE(sourceBias, "the fraction of photon packets distributed uniformly across secondary sources")
        ATTRIBUTE_MIN_VALUE(sourceBias, "[0")
        ATTRIBUTE_MAX_VALUE(sourceBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(sourceBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(sourceBias, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

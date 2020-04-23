/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYAOPTIONS_HPP
#define LYAOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The LyaOptions class simply offers a number of configuration options related to the treatment
    of Lyman-alpha line transfer, if this is enabled in the simulation. */
class LyaOptions : public SimulationItem
{
    /** The enumeration type indicating the supported Lyman-alpha acceleration schemes. */
    ENUM_DEF(LyaAccelerationScheme, None, Constant, Laursen, Smith)
        ENUM_VAL(LyaAccelerationScheme, None, "no acceleration")
        ENUM_VAL(LyaAccelerationScheme, Constant, "acceleration with a constant critical value")
        ENUM_VAL(LyaAccelerationScheme, Laursen, "the Laursen 2009 acceleration scheme")
        ENUM_VAL(LyaAccelerationScheme, Smith, "the Smith 2015 acceleration scheme")
    ENUM_END()

    ITEM_CONCRETE(LyaOptions, SimulationItem, "a set of options related to Lyman-alpha line transfer")

        PROPERTY_ENUM(lyaAccelerationScheme, LyaAccelerationScheme, "the Lyman-alpha line transfer acceleration scheme")
        ATTRIBUTE_DEFAULT_VALUE(lyaAccelerationScheme, "Constant")
        ATTRIBUTE_DISPLAYED_IF(lyaAccelerationScheme, "Level2")

        PROPERTY_DOUBLE(lyaAccelerationCriticalValue,
                        "the critical value to be employed by the 'constant' acceleration scheme")
        ATTRIBUTE_MIN_VALUE(lyaAccelerationCriticalValue, "]0")
        ATTRIBUTE_MAX_VALUE(lyaAccelerationCriticalValue, "100]")
        ATTRIBUTE_DEFAULT_VALUE(lyaAccelerationCriticalValue, "3")
        ATTRIBUTE_RELEVANT_IF(lyaAccelerationCriticalValue, "lyaAccelerationSchemeConstant")
        ATTRIBUTE_DISPLAYED_IF(lyaAccelerationCriticalValue, "Level2")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

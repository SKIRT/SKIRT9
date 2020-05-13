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
    ENUM_DEF(LyaAccelerationScheme, None, Constant, Variable)
        ENUM_VAL(LyaAccelerationScheme, None, "no acceleration")
        ENUM_VAL(LyaAccelerationScheme, Constant, "acceleration scheme with a constant critical value")
        ENUM_VAL(LyaAccelerationScheme, Variable, "acceleration scheme depending on local gas temperature and density")
    ENUM_END()

    ITEM_CONCRETE(LyaOptions, SimulationItem, "a set of options related to Lyman-alpha line transfer")

        PROPERTY_ENUM(lyaAccelerationScheme, LyaAccelerationScheme, "the Lyman-alpha line transfer acceleration scheme")
        ATTRIBUTE_DEFAULT_VALUE(lyaAccelerationScheme, "Variable")
        ATTRIBUTE_DISPLAYED_IF(lyaAccelerationScheme, "Level2")

        PROPERTY_DOUBLE(lyaAccelerationStrength, "the acceleration strength; higher is faster but less accurate")
        ATTRIBUTE_MIN_VALUE(lyaAccelerationStrength, "]0")
        ATTRIBUTE_MAX_VALUE(lyaAccelerationStrength, "10]")
        ATTRIBUTE_DEFAULT_VALUE(lyaAccelerationStrength, "1")
        ATTRIBUTE_RELEVANT_IF(lyaAccelerationStrength, "lyaAccelerationSchemeConstant|lyaAccelerationSchemeVariable")
        ATTRIBUTE_DISPLAYED_IF(lyaAccelerationStrength, "Level2")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

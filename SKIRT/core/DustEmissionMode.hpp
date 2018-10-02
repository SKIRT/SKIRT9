/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIONMODE_HPP
#define DUSTEMISSIONMODE_HPP

#include "WithMediumMode.hpp"

////////////////////////////////////////////////////////////////////

/** DustEmissionMode is a SimulationMode subclass indicating a simulation mode that includes
    secondary emission from dust, in addition to the effects of absorption and scattering. In this
    mode, the simulation keeps track of the radation field (to calculate the dust emission) and
    optionally performs iterations to self-consistently calculate the effects of dust
    self-absorption.

    This simulation mode is meaningful only for a continuous wavelength range that includes optical
    to far-infrared regimes. */
class DustEmissionMode : public WithMediumMode
{
    ITEM_CONCRETE(DustEmissionMode, WithMediumMode, "a simulation with secondary emission from dust")
        ATTRIBUTE_TYPE_ALLOWED_IF(DustEmissionMode, "Panchromatic")

    PROPERTY_DOUBLE(primaryPacketsMultiplier,
                    "the multiplier on the number of photon packets launched from primary sources")
        ATTRIBUTE_MIN_VALUE(primaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(primaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(primaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(primaryPacketsMultiplier, "Level3")

    PROPERTY_DOUBLE(secondaryPacketsMultiplier,
                    "the multiplier on the number of photon packets launched from secondary sources")
        ATTRIBUTE_MIN_VALUE(secondaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(secondaryPacketsMultiplier, "Level3")

    PROPERTY_INT(dummy, "dummy option for dust emission")
        ATTRIBUTE_DEFAULT_VALUE(dummy, "99")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

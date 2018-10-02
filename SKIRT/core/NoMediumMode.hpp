/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NOMEDIUMMODE_HPP
#define NOMEDIUMMODE_HPP

#include "SimulationMode.hpp"

////////////////////////////////////////////////////////////////////

/** NoMediumMode is a SimulationMode subclass indicating a simulation mode without any media, i.e. with
    primary sources only. It does not offer any additional options. */
class NoMediumMode : public SimulationMode
{
    ITEM_CONCRETE(NoMediumMode, SimulationMode, "a simulation without transfer medium (primary sources only)")
        ATTRIBUTE_TYPE_INSERT(NoMediumMode, "NoMedium")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

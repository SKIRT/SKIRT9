/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NOMEDIAMODE_HPP
#define NOMEDIAMODE_HPP

#include "SimulationMode.hpp"

////////////////////////////////////////////////////////////////////

/** NoMediaMode is a SimulationMode subclass indicating a simulation mode without media, i.e. with
    primary sources only. It offers just the basic configuration option defining the number of
    photon packets to be launched in the simulation. */
class NoMediaMode : public SimulationMode
{
    ITEM_CONCRETE(NoMediaMode, SimulationMode, "a simulation mode without media (primary sources only)")

    PROPERTY_DOUBLE(numPackets, "the default number of photon packets launched per simulation segment")
        ATTRIBUTE_MIN_VALUE(numPackets, "[0")
        ATTRIBUTE_MAX_VALUE(numPackets, "1e19]")
        ATTRIBUTE_DEFAULT_VALUE(numPackets, "1e6")

    PROPERTY_DOUBLE(primaryPacketsMultiplier,
                    "the multiplier on the number of photon packets launched from primary sources")
        ATTRIBUTE_MIN_VALUE(primaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(primaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(primaryPacketsMultiplier, "1")

    ITEM_END()

    /** \fn numPackets
        The number of photon packets is specified as a double-precision floating point number
        rather than as a 64-bit integer to avoid implementing yet another discoverable property
        type. As a side benefit, one can use exponential notation to specify a large number of
        photon packets. Also, note that a double can exactly represent all integers up to 9e15. The
        maximum number of photon packets is somewhat arbitrarily set to 1e19 because that number is
        close to the maximum number representable with a 64-bit unsigned integer. */
};

////////////////////////////////////////////////////////////////////

#endif

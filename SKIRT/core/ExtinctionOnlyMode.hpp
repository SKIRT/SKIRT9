/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EXTINCTIONONLYMODE_HPP
#define EXTINCTIONONLYMODE_HPP

#include "NoMediaMode.hpp"

////////////////////////////////////////////////////////////////////

/** ExtinctionOnlyMode is a SimulationMode subclass indicating a simulation mode that calculates
    the extinction of the primary radiation through the configured media, including the effects of
    absorption and scattering. In this mode, the simulation does not track the radation field nor
    the state of any media properties. This simulation mode is meaningful only for wavelengths at
    which secondary sources (radiation from the media) can be neglected, i.e. in the ultraviolet,
    optical and near-infrared. */
class ExtinctionOnlyMode : public NoMediaMode
{
    ITEM_CONCRETE(ExtinctionOnlyMode, NoMediaMode, "an extinction-only simulation mode (no secondary emission)")

    PROPERTY_DOUBLE(minWeightReduction, "the minimum weight reduction factor before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minWeightReduction, "[1e3")
        ATTRIBUTE_DEFAULT_VALUE(minWeightReduction, "1e4")

    PROPERTY_INT(minScattEvents, "the minimum number of forced scattering events before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minScattEvents, "0")
        ATTRIBUTE_MAX_VALUE(minScattEvents, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minScattEvents, "0")

    PROPERTY_DOUBLE(pathLengthBias,
                    "the fraction of path lengths sampled from a linear rather than an exponential distribution")
        ATTRIBUTE_MIN_VALUE(pathLengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(pathLengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(pathLengthBias, "0.5")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

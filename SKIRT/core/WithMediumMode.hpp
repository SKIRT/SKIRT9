/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WITHMEDIUMMODE_HPP
#define WITHMEDIUMMODE_HPP

#include "SimulationMode.hpp"

////////////////////////////////////////////////////////////////////

/** WithMediumMode is an abstract (i.e. intermediate) SimulationMode subclass indicating a
    simulation mode with a transfer medium. It offers options that are relevant as soon as there is
    a medium, e.g. related to the photon cycle and to sampling the medium density. */
class WithMediumMode : public SimulationMode
{
    ITEM_ABSTRACT(WithMediumMode, SimulationMode, "a simulation with a transfer medium")

    PROPERTY_DOUBLE(minWeightReduction, "the minimum weight reduction factor before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minWeightReduction, "[1e3")
        ATTRIBUTE_DEFAULT_VALUE(minWeightReduction, "1e4")
        ATTRIBUTE_DISPLAYED_IF(minWeightReduction, "Level3")

    PROPERTY_INT(minScattEvents, "the minimum number of forced scattering events before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minScattEvents, "0")
        ATTRIBUTE_MAX_VALUE(minScattEvents, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minScattEvents, "0")
        ATTRIBUTE_DISPLAYED_IF(minScattEvents, "Level3")

    PROPERTY_DOUBLE(pathLengthBias,
                    "the fraction of path lengths sampled from a stretched distribution")
        ATTRIBUTE_MIN_VALUE(pathLengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(pathLengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(pathLengthBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(pathLengthBias, "Level3")

    PROPERTY_INT(numDensitySamples, "the number of random density samples for determining spatial cell mass")
        ATTRIBUTE_MIN_VALUE(numDensitySamples, "10")
        ATTRIBUTE_MAX_VALUE(numDensitySamples, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numDensitySamples, "100")
        ATTRIBUTE_DISPLAYED_IF(numDensitySamples, "Level2")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

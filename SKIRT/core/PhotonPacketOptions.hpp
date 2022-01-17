/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTONPACKETOPTIONS_HPP
#define PHOTONPACKETOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The PhotonPacketOptions class simply offers a number of configuration options related to the
    Monte Carlo photon packet life cycle. These options are relevant as soon as there is a medium
    in the configuration.

    Two variations of the photon life cycle are implemented: with and without forced scattering.
    Forced scattering is the default behavior for all simulation modes except when Lyman-alpha line
    transfer is included (because the extra peel-offs for the many resonant scattering events tend
    to slow down the simulation). The implementation without forced scattering does \em not support
    storing the radiation field, which means it cannot be used when the simulation includes
    secondary emission or dynamic state iteration. */
class PhotonPacketOptions : public SimulationItem
{
    ITEM_CONCRETE(PhotonPacketOptions, SimulationItem, "a set of options related to the photon packet lifecycle")

        PROPERTY_BOOL(forceScattering, "use forced scattering to accelerate the photon life cycle")
        ATTRIBUTE_DEFAULT_VALUE(forceScattering, "Lya:false;true")
        ATTRIBUTE_RELEVANT_IF(forceScattering, "!(Emission|DynamicState)")
        ATTRIBUTE_DISPLAYED_IF(forceScattering, "Level3")
        ATTRIBUTE_INSERT(forceScattering, "Emission|DynamicState|forceScattering:ForceScattering")

        PROPERTY_DOUBLE(minWeightReduction, "the minimum weight reduction factor before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minWeightReduction, "[1e3")
        ATTRIBUTE_DEFAULT_VALUE(minWeightReduction, "1e4")
        ATTRIBUTE_RELEVANT_IF(minWeightReduction, "ForceScattering")
        ATTRIBUTE_DISPLAYED_IF(minWeightReduction, "Level3")

        PROPERTY_INT(minScattEvents,
                     "the minimum number of forced scattering events before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minScattEvents, "0")
        ATTRIBUTE_MAX_VALUE(minScattEvents, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minScattEvents, "0")
        ATTRIBUTE_RELEVANT_IF(minScattEvents, "ForceScattering")
        ATTRIBUTE_DISPLAYED_IF(minScattEvents, "Level3")

        PROPERTY_DOUBLE(pathLengthBias, "the fraction of path lengths sampled from a stretched distribution")
        ATTRIBUTE_MIN_VALUE(pathLengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(pathLengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(pathLengthBias, "Lya:0;0.5")
        ATTRIBUTE_RELEVANT_IF(pathLengthBias, "(ForceScattering)&(!Lya)")
        ATTRIBUTE_DISPLAYED_IF(pathLengthBias, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

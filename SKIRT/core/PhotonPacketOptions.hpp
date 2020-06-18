/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTONPACKETOPTIONS_HPP
#define PHOTONPACKETOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The PhotonPacketOptions class simply offers a number of configuration options related to the
    Monte Carlo photon packet lifecycle, such as when a photon packet should be terminated. These options
    are relevant as soon as there is a medium in the configuration. */
class PhotonPacketOptions : public SimulationItem
{
    ITEM_CONCRETE(PhotonPacketOptions, SimulationItem, "a set of options related to the photon packet lifecycle")

        PROPERTY_BOOL(forceScattering, "use forced scattering to accelerate the photon life cycle")
        ATTRIBUTE_DEFAULT_VALUE(forceScattering, "Lya:false;true")
        ATTRIBUTE_RELEVANT_IF(forceScattering, "ExtinctionOnly")
        ATTRIBUTE_DISPLAYED_IF(forceScattering, "Level3")
        ATTRIBUTE_INSERT(forceScattering, "forceScattering:ForceScattering")

        PROPERTY_DOUBLE(minWeightReduction, "the minimum weight reduction factor before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minWeightReduction, "[1e3")
        ATTRIBUTE_DEFAULT_VALUE(minWeightReduction, "1e4")
        ATTRIBUTE_RELEVANT_IF(minWeightReduction, "ForceScattering|Emission")
        ATTRIBUTE_DISPLAYED_IF(minWeightReduction, "Level3")

        PROPERTY_INT(minScattEvents,
                     "the minimum number of forced scattering events before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minScattEvents, "0")
        ATTRIBUTE_MAX_VALUE(minScattEvents, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minScattEvents, "0")
        ATTRIBUTE_RELEVANT_IF(minScattEvents, "ForceScattering|Emission")
        ATTRIBUTE_DISPLAYED_IF(minScattEvents, "Level3")

        PROPERTY_DOUBLE(pathLengthBias, "the fraction of path lengths sampled from a stretched distribution")
        ATTRIBUTE_MIN_VALUE(pathLengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(pathLengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(pathLengthBias, "0.5")
        ATTRIBUTE_RELEVANT_IF(pathLengthBias, "(ForceScattering|Emission)&(!Lya)")
        ATTRIBUTE_DISPLAYED_IF(pathLengthBias, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

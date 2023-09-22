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

    Several variations of the photon life cycle implementation can be configured:

    - With or without explicit absorption. Explicit absorption allows absorption cross sections to
    be negative, which can be the case for materials that exhibit stimulated emission. Note that
    scattering cross sections never can be negative. The default behavior is not to use explicit
    absorption, because the effects of this recent technique has not yet been well studied.

    - With or without forced scattering. Forced scattering tends to reduce noise for simulations
    with low to limited optical depths, such as for most dust models on galaxy-wide scales.
    Therefore, forced scattering is the default behavior except when Lyman-alpha line transfer is
    included, because the extra peel-offs for the many resonant scattering events tend to slow down
    the simulation. Furthermore, the implementation without forced scattering currently does \em
    not support storing the radiation field, which means it cannot be used when the simulation
    includes secondary emission or dynamic state iteration.

    The remaining three options serve to further configure the detailed behavior of the forced
    scattering photon cycle. The first two options determine when the photon life cycle will be
    terminated, i.e. after the weight of photon packet has decreased to some insignificant fraction
    of its original weight (luminosity) through interaction with the medium, with a given minimum
    number of scattering events if so desired.

    The last option configures the path length stretching mechanism. When determining the next
    interaction location of a photon packet with the medium, this technique samples in part from a
    distribution representing unphysically long path segments, correcting this deviation through a
    bias factor on the photon packet's weight. This allows a photon packet to more easily skip
    through high optical depth regions. The value configured for this option determines the
    fraction of path lengths sampled from the unphysical distribution. A value of zero disables the
    mechanism altogether.

    The path length stretching mechanism cannot be used for models where the photon packet
    wavelength may change during its life cycle, for example during scattering or as a result of
    bulk motions in the medium. This is so because the unphysically long distances between
    interactions would no longer correctly sample the effects on the wavelength. Also, path length
    stretching is currently \em not implemented for photon cycles without forced scattering. In all
    these cases, the path length stretching mechanism will automatically be disabled during setup.
    As a result, these simulations will lack the potential optimization brought by the path length
    technique. In particulatar, penetrating regions of high optical depth may require many
    scattering events with correspondingly longer running times. */
class PhotonPacketOptions : public SimulationItem
{
    ITEM_CONCRETE(PhotonPacketOptions, SimulationItem, "a set of options related to the photon packet lifecycle")

        PROPERTY_BOOL(explicitAbsorption, "use explicit absorption to allow negative absorption (stimulated emission)")
        ATTRIBUTE_DEFAULT_VALUE(explicitAbsorption, "false")
        ATTRIBUTE_DISPLAYED_IF(explicitAbsorption, "Level3")

        PROPERTY_BOOL(forceScattering, "use forced scattering to reduce noise")
        ATTRIBUTE_DEFAULT_VALUE(forceScattering, "Lya:false;true")
        ATTRIBUTE_RELEVANT_IF(forceScattering, "!(Emission|IteratePrimary)")
        ATTRIBUTE_DISPLAYED_IF(forceScattering, "Level3")
        ATTRIBUTE_INSERT(forceScattering, "Emission|IteratePrimary|forceScattering:ForceScattering")

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
        ATTRIBUTE_DEFAULT_VALUE(pathLengthBias, "(ForceScattering)&(!Lya):0.5;0")
        ATTRIBUTE_RELEVANT_IF(pathLengthBias, "(ForceScattering)&(!Lya)")
        ATTRIBUTE_DISPLAYED_IF(pathLengthBias, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

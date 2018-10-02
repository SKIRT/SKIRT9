/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATIONMODE_HPP
#define SIMULATIONMODE_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** SimulationMode is a simple class that allows the user to configure the fundamental operation
    mode of a SKIRT simulation and some related general options. The relevant options greatly
    depend on the mode of operation. For example, a simulation mode calculating the dust
    temperature self-consistently (taking into account dust self-absorption through iteration)
    would require more extensive configuration options than a mode restricted to calculating the
    extinction in the optical wavelength range.

    In practice, the SimulationMode class bundles configuration options related to the Monte Carlo
    photon lifecycle (e.g. the number of photon packets to be launched), the media state (e.g. the
    wavelength resolution of the radiation field grid), and the iterative simulation process (e.g.
    convergence criteria).

    The simulation mode should be included as one of the first properties of the top-level
    MonteCarloSimulation class, requiring the user to select a fundamental mode early in the
    configuration process. This in turn facilitates limiting the options offered during the
    remaining configuration process to those relevant for the selected mode. */
class SimulationMode : public SimulationItem
{
    /** The enumeration type indicating the treatment of media in the simulation. NoMedium
        indicates a simulation without transfer media, i.e. with primary sources only.
        ExtinctionOnly indicates a simulation that calculates the extinction of the primary
        radiation through the configured media, including the effects of scattering and absorption.
        In this mode, the simulation does not track the radation field nor the state of any media
        properties. This simulation mode is meaningful only for wavelengths at which secondary
        sources (radiation from the media) can be neglected, i.e. in the ultraviolet, optical and
        near-infrared. */
    ENUM_DEF(MediaTreatment, NoMedium, ExtinctionOnly)
    ENUM_VAL(MediaTreatment, NoMedium, "No transfer medium (primary sources only)")
    ENUM_VAL(MediaTreatment, ExtinctionOnly, "Include extinction (scattering and absorption); no secondary emission")
    ENUM_END()

    ITEM_CONCRETE(SimulationMode, SimulationItem, "a simulation mode")

    PROPERTY_ENUM(mediaTreatment, MediaTreatment, "the treatment of transfer media in the simulation")
        ATTRIBUTE_DEFAULT_VALUE(mediaTreatment, "ExtinctionOnly")
        ATTRIBUTE_INSERT(mediaTreatment, "mediaTreatmentNoMedium:NoMedium;"
                                         "mediaTreatmentExtinctionOnly:ExtinctionOnly")

    PROPERTY_DOUBLE(minWeightReduction, "the minimum weight reduction factor before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minWeightReduction, "[1e3")
        ATTRIBUTE_DEFAULT_VALUE(minWeightReduction, "1e4")
        ATTRIBUTE_RELEVANT_IF(minWeightReduction, "!NoMedium")
        ATTRIBUTE_DISPLAYED_IF(minWeightReduction, "Level3")

    PROPERTY_INT(minScattEvents, "the minimum number of forced scattering events before a photon packet is terminated")
        ATTRIBUTE_MIN_VALUE(minScattEvents, "0")
        ATTRIBUTE_MAX_VALUE(minScattEvents, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minScattEvents, "0")
        ATTRIBUTE_RELEVANT_IF(minScattEvents, "!NoMedium")
        ATTRIBUTE_DISPLAYED_IF(minScattEvents, "Level3")

    PROPERTY_DOUBLE(pathLengthBias,
                    "the fraction of path lengths sampled from a linear rather than an exponential distribution")
        ATTRIBUTE_MIN_VALUE(pathLengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(pathLengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(pathLengthBias, "0.5")
        ATTRIBUTE_RELEVANT_IF(pathLengthBias, "!NoMedium")
        ATTRIBUTE_DISPLAYED_IF(pathLengthBias, "Level3")

    PROPERTY_INT(numDensitySamples, "the number of random density samples for determining spatial cell mass")
        ATTRIBUTE_MIN_VALUE(numDensitySamples, "10")
        ATTRIBUTE_MAX_VALUE(numDensitySamples, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numDensitySamples, "100")
        ATTRIBUTE_RELEVANT_IF(numDensitySamples, "!NoMedium")
        ATTRIBUTE_DISPLAYED_IF(numDensitySamples, "Level2")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

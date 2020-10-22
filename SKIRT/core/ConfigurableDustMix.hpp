/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURABLEDUSTMIX_HPP
#define CONFIGURABLEDUSTMIX_HPP

#include "GrainPopulation.hpp"
#include "MultiGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The ConfigurableDustMix class represents a fully user-configurable dust mix described by one or
    more dust grain populations. Specifically, the class can be configured with a list of
    GrainPopulation instances, each of which represents a particular dust grain population with
    configurable grain composition, grain size distribution, and size bin discretization. During
    setup, this list of grain populations is simply passed to the base class. */
class ConfigurableDustMix : public MultiGrainDustMix
{
    /** The enumeration type indicating the scattering mode. */
    ENUM_DEF(ScatteringType, HenyeyGreenstein, MaterialPhaseFunction, SphericalPolarization, SpheroidalPolarization)
        ENUM_VAL(ScatteringType, HenyeyGreenstein, "use the Henyey-Greenstein phase function (unpolarized)")
        ENUM_VAL(ScatteringType, MaterialPhaseFunction,
                 "use the phase function derived from actual material properties (unpolarized)")
        ENUM_VAL(ScatteringType, SphericalPolarization, "support polarization through scattering by spherical grains")
        ENUM_VAL(ScatteringType, SpheroidalPolarization,
                 "support polarization through scattering, absorption and emission by spheroidal grains")
    ENUM_END()

    ITEM_CONCRETE(ConfigurableDustMix, MultiGrainDustMix, "a configurable dust mix with one or more grain populations")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ConfigurableDustMix, "Level2")

        PROPERTY_ENUM(scatteringType, ScatteringType, "the type of scattering to be implemented")
        ATTRIBUTE_DEFAULT_VALUE(scatteringType, "HenyeyGreenstein")
        ATTRIBUTE_INSERT(scatteringType, "scatteringTypeSpheroidalPolarization:Spheroidal")

        PROPERTY_ITEM_LIST(populations, GrainPopulation, "the grain populations")
        ATTRIBUTE_DEFAULT_VALUE(populations, "GrainPopulation")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the configured grain populations to the dust mix. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns the scattering mode supported by this material mix as configured by
        the user through the scatteringType property. */
    ScatteringMode scatteringMode() const override;

    /** This function returns a flag indicating whether the material mix supports polarization
        during scattering events or not. For this dust mix, the function returns true if the user
        configured the SphericalPolarization or the SpheroidalPolarization scattering type, and
        false otherwise. */
    bool hasPolarizedScattering() const override;

    /** This function returns a flag indicating whether the secondary emission for this material
        mix may be polarized and anisotropic. For this dust mix, the function returns true if the
        user configured the SpheroidalPolarization scattering type, and false otherwise. */
    bool hasPolarizedEmission() const override;
};

////////////////////////////////////////////////////////////////////

#endif

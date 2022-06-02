/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EMITTINGGASMIX_HPP
#define EMITTINGGASMIX_HPP

#include "MaterialMix.hpp"
#include "SourceWavelengthRangeInterface.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** The EmittingGasMix is an abstract base class describing gas material mixes that support
    secondary line or continuum emission. It provides generic configuration options related to
    secondary gas emission. */
class EmittingGasMix : public MaterialMix, public SourceWavelengthRangeInterface
{
    ITEM_ABSTRACT(EmittingGasMix, MaterialMix, "an emitting gas mix")
        ATTRIBUTE_TYPE_INSERT(EmittingGasMix, "GasMix")

        PROPERTY_DOUBLE(sourceWeight, "the weight of this secondary source for the number of photon packets launched")
        ATTRIBUTE_MIN_VALUE(sourceWeight, "]0")
        ATTRIBUTE_MAX_VALUE(sourceWeight, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(sourceWeight, "1")
        ATTRIBUTE_DISPLAYED_IF(sourceWeight, "Level3")

        PROPERTY_DOUBLE(wavelengthBias,
                        "the fraction of secondary photon packet wavelengths sampled from a bias distribution")
        ATTRIBUTE_MIN_VALUE(wavelengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(wavelengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBias, "Level3")

        PROPERTY_ITEM(wavelengthBiasDistribution, WavelengthDistribution,
                      "the bias distribution for sampling secondary photon packet wavelengths")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBiasDistribution, "LogWavelengthDistribution")
        ATTRIBUTE_RELEVANT_IF(wavelengthBiasDistribution, "wavelengthBias")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBiasDistribution, "Level3")

    ITEM_END()

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

    //======== Wavelength range interface =======

public:
    /** This function returns the wavelength range of the gas emission wavelength grid,
        implementing the SourceWavelengthRangeInterface interface so that the wavelength bias
        distribution configured here can locate its associated "source" wavelength range.

        The function determines and returns the union of wavelength ranges published by the
        subclass for line and continuum emission. */
    Range wavelengthRange() const override;
};

////////////////////////////////////////////////////////////////////

#endif

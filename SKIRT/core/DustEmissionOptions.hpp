/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIONOPTIONS_HPP
#define DUSTEMISSIONOPTIONS_HPP

#include "SimulationItem.hpp"
#include "DustEmissivity.hpp"
#include "WavelengthDistribution.hpp"
#include "DisjointWavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** The DustEmissionOptions class simply offers a number of configuration options related to
    secondary emission from dust. In a mode where dust emission is enabled, the simulation keeps
    track of the radation field and it also needs a wavelength grid on which to calculate the dust
    emission spectrum. */
class DustEmissionOptions : public SimulationItem
{
    ITEM_CONCRETE(DustEmissionOptions, SimulationItem, "a set of options related to secondary emission from dust")

    PROPERTY_ITEM(dustEmissivity, DustEmissivity, "the dust emissivity calculator")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissivity, "EquilibriumDustEmissivity")

//    PROPERTY_ITEM(cellLibrary, CellLibrary, "the library mechanism for combining spatial cells")
//        ATTRIBUTE_DEFAULT_VALUE(cellLibrary, "AllCellsLibrary")

    PROPERTY_ITEM(radiationFieldWLG, DisjointWavelengthGrid, "the wavelength grid for storing the radiation field")
        ATTRIBUTE_DEFAULT_VALUE(radiationFieldWLG, "LogWavelengthGrid")

    PROPERTY_ITEM(dustEmissionWLG, DisjointWavelengthGrid, "the wavelength grid for calculating the dust emission spectrum")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissionWLG, "LogWavelengthGrid")

    PROPERTY_DOUBLE(secondaryPacketsMultiplier,
                    "the multiplier on the number of photon packets launched for secondary emission from dust")
        ATTRIBUTE_MIN_VALUE(secondaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(secondaryPacketsMultiplier, "Level3")

    PROPERTY_DOUBLE(spatialBias, "the fraction of secondary photon packets distributed uniformly across spatial cells")
        ATTRIBUTE_MIN_VALUE(spatialBias, "[0")
        ATTRIBUTE_MAX_VALUE(spatialBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(spatialBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(spatialBias, "Level3")

    PROPERTY_DOUBLE(wavelengthBias,
                    "the fraction of secondary photon packet wavelengths sampled from a bias distribution")
        ATTRIBUTE_MIN_VALUE(wavelengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(wavelengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBias, "0.5")
        ATTRIBUTE_RELEVANT_IF(wavelengthBias, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBias, "Level3")

    PROPERTY_ITEM(wavelengthBiasDistribution, WavelengthDistribution,
                  "the bias distribution for sampling secondary photon packet wavelengths")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBiasDistribution, "LogWavelengthDistribution")
        ATTRIBUTE_RELEVANT_IF(wavelengthBiasDistribution, "wavelengthBias")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBiasDistribution, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

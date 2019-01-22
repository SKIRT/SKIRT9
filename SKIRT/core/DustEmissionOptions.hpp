/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIONOPTIONS_HPP
#define DUSTEMISSIONOPTIONS_HPP

#include "SimulationItem.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "SpatialCellLibrary.hpp"
#include "WavelengthDistribution.hpp"
#include "WavelengthRangeInterface.hpp"

////////////////////////////////////////////////////////////////////

/** The DustEmissionOptions class simply offers a number of configuration options related to
    secondary emission from dust. In a mode where dust emission is enabled, the simulation keeps
    track of the radation field and it also needs a wavelength grid on which to calculate the dust
    emission spectrum. */
class DustEmissionOptions : public SimulationItem, public WavelengthRangeInterface
{
    /** The enumeration type indicating the method used for dust emission calculations. */
    ENUM_DEF(EmissionType, Equilibrium, Stochastic)
    ENUM_VAL(EmissionType, Equilibrium, "assume local thermal equilibrium (LTE)")
    ENUM_VAL(EmissionType, Stochastic, "handle stochastically heated grains (non-LTE)")
    ENUM_END()

    ITEM_CONCRETE(DustEmissionOptions, SimulationItem, "a set of options related to secondary emission from dust")

    PROPERTY_ENUM(dustEmissionType, EmissionType, "the method used for dust emission calculations")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissionType, "Equilibrium")
        ATTRIBUTE_INSERT(dustEmissionType, "dustEmissionTypeStochastic:StochasticDustEmission")
        ATTRIBUTE_DISPLAYED_IF(dustEmissionType, "Level2")

    PROPERTY_ITEM(cellLibrary, SpatialCellLibrary, "the spatial cell grouping scheme for calculating dust emission")
        ATTRIBUTE_DEFAULT_VALUE(cellLibrary, "AllCellsLibrary")
        ATTRIBUTE_REQUIRED_IF(cellLibrary, "false")
        ATTRIBUTE_DISPLAYED_IF(cellLibrary, "Level2")

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

    //======================== Other Functions =======================

public:
    /** This function returns the wavelength range of the dust emission wavelength grid,
        implementing the WavelengthRangeInterface interface so that the wavelength bias
        distribution configured here can locate its associated "source" wavelength range. */
    Range wavelengthRange() const override { return _dustEmissionWLG->wavelengthRange(); }
};

////////////////////////////////////////////////////////////////////

#endif

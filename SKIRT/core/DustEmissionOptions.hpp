/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIONOPTIONS_HPP
#define DUSTEMISSIONOPTIONS_HPP

#include "DisjointWavelengthGrid.hpp"
#include "SimulationItem.hpp"
#include "SourceWavelengthRangeInterface.hpp"
#include "SpatialCellLibrary.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** The DustEmissionOptions class simply offers a number of configuration options related to
    thermal emission from dust. In a mode where dust emission is enabled, the simulation also
    needs a wavelength grid on which to calculate the dust emission spectrum. */
class DustEmissionOptions : public SimulationItem, public SourceWavelengthRangeInterface
{
    /** The enumeration type indicating the method used for dust emission calculations. */
    ENUM_DEF(EmissionType, Equilibrium, Stochastic)
        ENUM_VAL(EmissionType, Equilibrium, "assume local thermal equilibrium (LTE)")
        ENUM_VAL(EmissionType, Stochastic, "handle stochastically heated grains (non-LTE)")
    ENUM_END()

    ITEM_CONCRETE(DustEmissionOptions, SimulationItem, "a set of options related to thermal emission from dust")

        PROPERTY_ENUM(dustEmissionType, EmissionType, "the method used for dust emission calculations")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissionType, "Equilibrium")
        ATTRIBUTE_INSERT(dustEmissionType, "dustEmissionTypeStochastic:StochasticDustEmission")
        ATTRIBUTE_DISPLAYED_IF(dustEmissionType, "Level2")

        PROPERTY_BOOL(includeHeatingByCMB, "add the cosmic microwave background (CMB) as a dust heating source term")
        ATTRIBUTE_DEFAULT_VALUE(includeHeatingByCMB, "false")
        ATTRIBUTE_DISPLAYED_IF(includeHeatingByCMB, "NonZeroRedshift")

        PROPERTY_ITEM(cellLibrary, SpatialCellLibrary, "the spatial cell grouping scheme for calculating dust emission")
        ATTRIBUTE_DEFAULT_VALUE(cellLibrary, "AllCellsLibrary")
        ATTRIBUTE_REQUIRED_IF(cellLibrary, "false")
        ATTRIBUTE_DISPLAYED_IF(cellLibrary, "Level2")

        PROPERTY_ITEM(dustEmissionWLG, DisjointWavelengthGrid,
                      "the wavelength grid for calculating the dust emission spectrum")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissionWLG, "LogWavelengthGrid")

        PROPERTY_DOUBLE(maxFractionOfPrimary, "convergence is reached when the total absorbed dust luminosity "
                                              "is less than this fraction of the total absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrimary, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrimary, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrimary, "0.01")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrimary, "IterateSecondary")

        PROPERTY_DOUBLE(maxFractionOfPrevious,
                        "convergence is reached when the total absorbed dust luminosity "
                        "has changed by less than this fraction compared to the previous iteration")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrevious, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrevious, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrevious, "0.03")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrevious, "IterateSecondary")

        PROPERTY_DOUBLE(sourceWeight, "the weight of dust emission for the number of photon packets launched")
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

    //======================== Other Functions =======================

public:
    /** This function returns the wavelength range of the dust emission wavelength grid,
        implementing the SourceWavelengthRangeInterface interface so that the wavelength bias
        distribution configured here can locate its associated "source" wavelength range. */
    Range wavelengthRange() const override { return _dustEmissionWLG->wavelengthRange(); }
};

////////////////////////////////////////////////////////////////////

#endif

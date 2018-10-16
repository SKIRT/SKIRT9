/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIONMODE_HPP
#define DUSTEMISSIONMODE_HPP

#include "WithMediumMode.hpp"
#include "DustEmissivity.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** DustEmissionMode is a SimulationMode subclass indicating a simulation mode that includes
    secondary emission from dust, in addition to the effects of absorption and scattering. In this
    mode, the simulation keeps track of the radation field (to calculate the dust emission) and
    optionally performs iterations to self-consistently calculate the effects of dust
    self-absorption.

    This simulation mode is meaningful only for a continuous wavelength range that includes optical
    to far-infrared regimes. */
class DustEmissionMode : public WithMediumMode
{
    ITEM_CONCRETE(DustEmissionMode, WithMediumMode, "a simulation with secondary emission from dust")
        ATTRIBUTE_TYPE_ALLOWED_IF(DustEmissionMode, "Panchromatic")
        ATTRIBUTE_TYPE_INSERT(DustEmissionMode, "DustEmission,Emission,RadiationField")

    PROPERTY_ITEM(dustEmissivity, DustEmissivity, "the dust emissivity calculator")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissivity, "EquilibriumDustEmissivity")

//    PROPERTY_ITEM(cellLibrary, CellLibrary, "the library mechanism for combining spatial cells")
//        ATTRIBUTE_DEFAULT_VALUE(cellLibrary, "AllCellsLibrary")

    PROPERTY_ITEM(radiationFieldWLG, WavelengthGrid, "the wavelength grid for storing the radiation field")
        ATTRIBUTE_DEFAULT_VALUE(radiationFieldWLG, "LogWavelengthGrid")

    PROPERTY_ITEM(emissionSpectrumWLG, WavelengthGrid, "the wavelength grid for calculating the dust emission spectrum")
        ATTRIBUTE_DEFAULT_VALUE(emissionSpectrumWLG, "LogWavelengthGrid")

    PROPERTY_DOUBLE(emissionBias, "the fraction of secondary photon packets distributed uniformly across spatial cells")
        ATTRIBUTE_MIN_VALUE(emissionBias, "[0")
        ATTRIBUTE_MAX_VALUE(emissionBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(emissionBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(emissionBias, "Level3")

    PROPERTY_BOOL(iterateSelfAbsorption, "self-consistently calculate dust self-absorption through iteration")
        ATTRIBUTE_DEFAULT_VALUE(iterateSelfAbsorption, "false")
        ATTRIBUTE_DISPLAYED_IF(iterateSelfAbsorption, "Level2")

    PROPERTY_INT(minIterations, "the minimum number of dust self-absorption iterations")
        ATTRIBUTE_MIN_VALUE(minIterations, "1")
        ATTRIBUTE_MAX_VALUE(minIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minIterations, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(minIterations, "Level3")

    PROPERTY_INT(maxIterations, "the maximum number of dust self-absorption iterations")
        ATTRIBUTE_MIN_VALUE(maxIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxIterations, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(maxIterations, "Level3")

    PROPERTY_DOUBLE(maxFractionOfPrimary, "convergence is reached when the total absorbed dust luminosity "
                                          "is less than this fraction of the total absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrimary, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrimary, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrimary, "0.01")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrimary, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrimary, "Level3")

    PROPERTY_DOUBLE(maxFractionOfPrevious, "convergence is reached when the total absorbed dust luminosity "
                                           "has changed by less than this fraction compared to the previous iteration")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrevious, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrevious, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrevious, "0.03")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrevious, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrevious, "Level3")

    PROPERTY_DOUBLE(primaryPacketsMultiplier,
                    "the multiplier on the number of photon packets launched from primary sources")
        ATTRIBUTE_MIN_VALUE(primaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(primaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(primaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(primaryPacketsMultiplier, "Level3")

    PROPERTY_DOUBLE(iterationPacketsMultiplier,
                    "the multiplier on the number of photon packets launched for each self-absorption iteration")
        ATTRIBUTE_MIN_VALUE(iterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(iterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(iterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(iterationPacketsMultiplier, "iterateSelfAbsorption")
        ATTRIBUTE_DISPLAYED_IF(iterationPacketsMultiplier, "Level3")

    PROPERTY_DOUBLE(secondaryPacketsMultiplier,
                    "the multiplier on the number of photon packets launched for secondary emission")
        ATTRIBUTE_MIN_VALUE(secondaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(secondaryPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

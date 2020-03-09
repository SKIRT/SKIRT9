/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASEMISSIONOPTIONS_HPP
#define GASEMISSIONOPTIONS_HPP

#include "DisjointWavelengthGrid.hpp"
#include "ItemInfo.hpp"
#include "SimulationItem.hpp"
#include "SourceWavelengthRangeInterface.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** Documentation */

class GasEmissionOptions : public SimulationItem, public SourceWavelengthRangeInterface
{
    ITEM_CONCRETE(GasEmissionOptions, SimulationItem, "a set of options related to secondary emission from gas")

        PROPERTY_ITEM(gasEmissionWLG, DisjointWavelengthGrid,
                      "the wavelength grid for calculating the gas emission spectrum")
        ATTRIBUTE_DEFAULT_VALUE(gasEmissionWLG, "LogWavelengthGrid")

        // there are a couple of features in DustEmissionOptions which also apply here, but should
        // maybe be moved into a common "SecondarayEmissionOptions" thing:

        // packet multiplier

        // store radiation field

        // spatial bias

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
    /** This function returns the wavelength range of the gas emission wavelength grid,
        implementing the SourceWavelengthRangeInterface interface so that the wavelength bias
        distribution configured here can locate its associated "source" wavelength range. */
    Range wavelengthRange() const override { return _gasEmissionWLG->wavelengthRange(); }
};

////////////////////////////////////////////////////////////////////

#endif

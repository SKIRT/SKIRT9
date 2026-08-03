/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECIALTYWAVELENGTHGRIDPROBE_HPP
#define SPECIALTYWAVELENGTHGRIDPROBE_HPP

#include "MaterialWavelengthRangeInterface.hpp"
#include "SpecialtyProbe.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** SpecialtyWavelengthGridProbe is a base class for probes that require a configurable wavelength
    grid property. It implements the MaterialWavelengthRangeInterface to indicate that
    wavelength-dependent material properties may be required for the configured wavelength grid. */
class SpecialtyWavelengthGridProbe : public SpecialtyProbe, public MaterialWavelengthRangeInterface
{
    ITEM_ABSTRACT(SpecialtyWavelengthGridProbe, SpecialtyProbe, "a specialty wavelength grid probe")

        PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for this probe")
        ATTRIBUTE_RELEVANT_IF(wavelengthGrid, "Panchromatic")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultInstrumentWavelengthGrid")
        ATTRIBUTE_DISPLAYED_IF(wavelengthGrid, "Level2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns a wavelength range corresponding to the user-configured wavelength
        grid for this probe, if any, indicating that wavelength-dependent material properties may be required for
        this wavelength range. */
    Range wavelengthRange() const override;

    /** This function returns a pointer to the user-configured wavelength grid for this probe, if
        any, indicating that wavelength-dependent material properties may be required for these
        wavelengths. */
    WavelengthGrid* materialWavelengthGrid() const override;
};

////////////////////////////////////////////////////////////////////

#endif

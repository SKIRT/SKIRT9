/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSTRACTWAVELENGTHGRIDPROBE_HPP
#define ABSTRACTWAVELENGTHGRIDPROBE_HPP

#include "Probe.hpp"
#include "MaterialWavelengthRangeInterface.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** AbstractWavelengthGridProbe is a base class for probes that require a configurable wavelength
    grid property. It implements the MaterialWavelengthRangeInterface to indicate that
    wavelength-dependent material properties may be required for the configured wavelength grid. */
class AbstractWavelengthGridProbe : public Probe, public MaterialWavelengthRangeInterface
{
    ITEM_ABSTRACT(AbstractWavelengthGridProbe, Probe, "a probe requiring a wavelength grid")

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

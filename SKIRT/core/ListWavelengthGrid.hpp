/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTWAVELENGTHGRID_HPP
#define LISTWAVELENGTHGRID_HPP

#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** ListWavelengthGrid is a subclass of the WavelengthGrid class that allows the list of
    representative wavelengths in the grid to be specified as a property in the configuration file.
    It is intended for use in cases where there are just a few wavelengths, but nothing keeps the
    user from specifying a long list. */
class ListWavelengthGrid : public WavelengthGrid
{
    ITEM_CONCRETE(ListWavelengthGrid, WavelengthGrid, "a list of one or more wavelengths")

    PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the representative wavelengths specified by the
        \em wavelengths property. */
    void setupSelfBefore() override;
};

//////////////////////////////////////////////////////////////////////

#endif

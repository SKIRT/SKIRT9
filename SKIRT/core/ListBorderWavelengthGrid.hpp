/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTBORDERWAVELENGTHGRID_HPP
#define LISTBORDERWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** ListBorderWavelengthGrid is a subclass of the DisjointWavelengthGrid class that allows an
    arbitrary wavelength grid to be fully specified inside the configuration file (i.e. without
    referring to an input file). It is intended for use in cases where there are just a few
    wavelengths, but nothing keeps the user from specifying a long list. The order of the specified
    wavelengths is not important; they will be sorted anyway.

    The wavelengths specified in the list represent the borders of the wavelength grid bins; there
    must be at least two borders in the list. The characteristic wavelength for each bin is derived
    automatically in linear or logarithmic space depending on the value of the \em log flag. For
    more details, refer to the DisjointWavelengthGrid::setWavelengthBorder() function. */
class ListBorderWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(ListBorderWavelengthGrid, DisjointWavelengthGrid,
                  "a wavelength grid configured as a list of bin borders")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ListBorderWavelengthGrid, "Level2")

        PROPERTY_DOUBLE_LIST(wavelengths, "the wavelength bin borders")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_BOOL(log, "logarithmic scale")
        ATTRIBUTE_DEFAULT_VALUE(log, "true")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the bin borders specified by the \em wavelengths
        property. */
    void setupSelfBefore() override;
};

//////////////////////////////////////////////////////////////////////

#endif

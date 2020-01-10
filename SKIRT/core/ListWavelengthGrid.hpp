/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTWAVELENGTHGRID_HPP
#define LISTWAVELENGTHGRID_HPP

#include "DisjointWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** ListWavelengthGrid is a subclass of the DisjointWavelengthGrid class that allows an arbitrary
    wavelength grid to be fully specified inside the configuration file (i.e. without referring to
    an input file). It is intended for use in cases where there are just a few wavelengths, but
    nothing keeps the user from specifying a long list. The order of the specified wavelengths is
    not important; they will be sorted anyway.

    The class constructs either a consecutive range of adjacent wavelength bins (when \em
    relativeHalfWidth is zero, the default value) or a set of distinct nonadjacent wavelength bins
    (with the relative half bin width given by a nonzero value of \em relativeHalfWidth). Refer to
    the DisjointWavelengthGrid class for more details. */
class ListWavelengthGrid : public DisjointWavelengthGrid
{
    ITEM_CONCRETE(ListWavelengthGrid, DisjointWavelengthGrid, "a list of one or more wavelengths")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ListWavelengthGrid, "Level2")

        PROPERTY_DOUBLE_LIST(wavelengths, "the characteristic wavelength for each bin")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 Angstrom")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_DOUBLE(relativeHalfWidth, "the relative half width for discrete bins, or zero for a consecutive range")
        ATTRIBUTE_MIN_VALUE(relativeHalfWidth, "[0")
        ATTRIBUTE_MAX_VALUE(relativeHalfWidth, "1[")
        ATTRIBUTE_DEFAULT_VALUE(relativeHalfWidth, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function sets the wavelength grid to the characteristic wavelengths specified by the
        \em wavelengths property. */
    void setupSelfBefore() override;
};

//////////////////////////////////////////////////////////////////////

#endif

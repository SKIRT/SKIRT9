/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILESED_HPP
#define FILESED_HPP

#include "TabulatedSED.hpp"

////////////////////////////////////////////////////////////////////

/** A FileSED object represents a spectral energy distribution that is loaded from an input file.
    The floating point numbers in the first two columns of the text file specify respectively the
    wavelength and the corresponding specific luminosity per unit of wavelength. Any additional
    columns in the file are ignored. The wavelengths must be listed in increasing or decreasing
    order. The luminosity outside the range indicated by the first and the last wavelength in the
    file is considered to be zero.

    The wavelengths are by default given in micron (the units can be overridden by column header
    info in the file). The specific luminosity values are by default given in per-wavelength units.
    Again, this can be overridden by column header info in the file to use per-frequency or neutral
    units. Other than this, the scaling of the values is arbitrary because the %SED will be
    normalized after being loaded. However, the input procedure still insists on knowing the
    precise units. */
class FileSED : public TabulatedSED
{
    ITEM_CONCRETE(FileSED, TabulatedSED, "a spectral energy distribution loaded from a text file")

        PROPERTY_STRING(filename, "the name of the file with the spectral energy distribution")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file into the specified arrays. */
    void getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

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
    columns in the file are ignored. The wavelengths are by default given in micron (the units can
    be overridden by column header info in the file) and must listed be in increasing order. The
    specific luminosity values can be given in arbitrary units because the %SED will be normalized
    after being loaded. The luminosity outside the range indicated by the first and the last
    wavelength in the file is considered to be zero. */
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

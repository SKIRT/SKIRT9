/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEBAND_HPP
#define FILEBAND_HPP

#include "TabulatedBand.hpp"

////////////////////////////////////////////////////////////////////

/** A FileBand object represents a wavelength band with a transmission curve that is loaded from an
    input file. The floating point numbers in the first two columns of the text file specify
    respectively the wavelength and the corresponding transmission value. Any additional columns in
    the file are ignored. The wavelengths must be listed in increasing or decreasing order. The
    transmission outside the range indicated by the first and the last wavelength in the list is
    considered to be zero.

    The wavelengths are by default given in micron (the units can be overridden by column header
    info in the file). The transmission values are dimensionless. The scaling of these values is
    arbitrary because the transmission curve will be normalized after being loaded. Refer to the
    description of the Band class for more information. */
class FileBand : public TabulatedBand
{
    ITEM_CONCRETE(FileBand, TabulatedBand, "a wavelength band (transmission curve) loaded from a text file")

        PROPERTY_STRING(filename, "the name of the file defining the transmission curve")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the wavelengths and the corresponding transmission values from the
        user-configured text file and returns them into the result arrays. */
    void getWavelengthsAndTransmissions(Array& lambdav, Array& transv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

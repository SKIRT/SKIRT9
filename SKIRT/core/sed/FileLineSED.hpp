/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILELINESED_HPP
#define FILELINESED_HPP

#include "LineSED.hpp"

////////////////////////////////////////////////////////////////////

/** A FileLineSED object represents a spectral energy distribution that consist of one or more
    discrete emission lines with zero width and is loaded from an input file. The floating point
    numbers in the first two columns of the text file specify respectively the line wavelengths and
    the corresponding luminosities. Any additional columns in the file are ignored.

    The wavelengths are by default given in micron (the units can be overridden by column header
    info in the file) and may be listed in any order. The corresponding luminosity values are by
    default given in solar units, which can be overridden by column header info in the file.
    However, although the input procedure insists on knowing the precise units, the scaling of the
    values is arbitrary because the %SED will be normalized after being loaded.

    Because the specific luminosity of the spectrum is undefined, there are important restrictions
    on the use of this %SED. See the description of the LineSED class for more information. */
class FileLineSED : public LineSED
{
    ITEM_CONCRETE(FileLineSED, LineSED, "a discrete line SED loaded from a text file")

        PROPERTY_STRING(filename, "the name of the file with the line definitions")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file into the specified arrays. */
    void getWavelengthsAndLuminosities(Array& lambdav, Array& Lv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

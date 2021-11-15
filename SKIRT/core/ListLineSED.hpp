/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTLINESED_HPP
#define LISTLINESED_HPP

#include "LineSED.hpp"

////////////////////////////////////////////////////////////////////

/** A ListLineSED object represents a spectral energy distribution that consist of one or more
    discrete emission lines with zero width and that is fully specified inside the configuration
    file (i.e. without referring to an input file). It is intended for use in cases where there are
    just a few wavelength/luminosity pairs, but nothing keeps the user from specifying a long list.

    The wavelengths are by default given in micron and may be listed be in any order. The
    corresponding luminosity values are by default given in solar luminosity units, but the scaling
    of the values is arbitrary because the %SED will be normalized after being loaded.

    Because the specific luminosity of the spectrum is undefined, there are important restrictions
    on the use of this %SED. See the description of the LineSED class for more information. */
class ListLineSED : public LineSED
{
    ITEM_CONCRETE(ListLineSED, LineSED, "a discrete line SED specified inside the configuration file")

        PROPERTY_DOUBLE_LIST(wavelengths, "the line wavelengths")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_DOUBLE_LIST(luminosities, "the line luminosities")
        ATTRIBUTE_QUANTITY(luminosities, "bolluminosity")
        ATTRIBUTE_MIN_VALUE(luminosities, "[0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the number of configured wavelengths and luminosities match and
        loads the data into the specified arrays. */
    void getWavelengthsAndLuminosities(Array& lambdav, Array& Lv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

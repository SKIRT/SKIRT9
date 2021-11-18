/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTBAND_HPP
#define LISTBAND_HPP

#include "TabulatedBand.hpp"

////////////////////////////////////////////////////////////////////

/** A ListBand object represents a wavelength band with a transmission curve that is fully
    specified inside the configuration file (i.e. without referring to an input file). It is
    intended for use in cases where there are just a few wavelength/transmission pairs, but nothing
    keeps the user from specifying a long list. The wavelengths must be listed in increasing or
    decreasing order, and a dimensionless transmission value must be specified corresponding to
    each wavelength. The transmission outside the range indicated by the first and the last
    wavelength in the list is considered to be zero. The scaling of the transmission values is
    arbitrary because the transmission curve will be normalized after being loaded. Refer to the
    description of the Band class for more information. */
class ListBand : public TabulatedBand
{
    ITEM_CONCRETE(ListBand, TabulatedBand,
                  "a wavelength band (transmission curve) specified inside the configuration file")

        PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths at which to specify the transmission")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

        PROPERTY_DOUBLE_LIST(transmissionValues, "the transmission at each of the given wavelengths")
        ATTRIBUTE_MIN_VALUE(transmissionValues, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function copies the user-configured wavelengths and transmission values into the
        result arrays. */
    void getWavelengthsAndTransmissions(Array& lambdav, Array& transv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LISTBAND_HPP
#define LISTBAND_HPP

#include "Band.hpp"

////////////////////////////////////////////////////////////////////

/** A ListBand object represents a wavelength band with a transmission curve that is fully
    specified inside the configuration file (i.e. without referring to an input file). It is
    intended for use in cases where there are just a few wavelength/transmission pairs, but nothing
    keeps the user from specifying a long list. The transmission outside the range indicated by the
    first and the last wavelength in the list is considered to be zero.

    The wavelengths must listed be in increasing order, and a transmission value must be specified
    corresponding to each wavelength. The scaling of the transmission values is arbitrary the
    transmission curve will be normalized after being loaded. Refer to the description of the Band
    class for more information. */
class ListBand : public Band
{
    ITEM_CONCRETE(ListBand, Band, "a wavelength band (transmission curve) specified inside the configuration file")

    PROPERTY_DOUBLE_LIST(wavelengths, "the wavelengths at which to specify the transmission")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")

    PROPERTY_DOUBLE_LIST(transmissionValues, "the transmission at each of the given wavelengths")
        ATTRIBUTE_MIN_VALUE(transmissionValues, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function normalizes the transmission values configured by the user as described in the
        Band class header. */
    void setupSelfBefore() override;

    //============= Functions required by base class =============

protected:
    /** This function returns the number of elements in the wavelength and transmission data arrays
        held by a subclass. */
    size_t dataSize() const override;

    /** This function returns a pointer to the first wavelength in the corresponding array held by
        a subclass. The number of elements in this array can be obtained through the dataSize()
        function. */
    const double* wavelengthData() const override;

    /** This function returns a pointer to the first transmission value in the corresponding array
        held by a subclass. The number of elements in this array can be obtained through the
        dataSize() function. The transmission values are normalized as described in the Band class
        header. */
    const double* transmissionData() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _transv;  // normalized transmission values
};

////////////////////////////////////////////////////////////////////

#endif

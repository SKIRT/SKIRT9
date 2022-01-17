/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABULATEDBAND_HPP
#define TABULATEDBAND_HPP

#include "Band.hpp"

////////////////////////////////////////////////////////////////////

/** TabulatedBand is an abstract class for representing a wavelength band that is tabulated by the
    user in the form of wavelength/transmission pairs. The wavelengths must be listed in increasing
    or decreasing order. The transmission outside the range indicated by the first and the last
    wavelength in the list is considered to be zero. Refer to the description of the Band class for
    more information.

    The subclass must load the tabulated data. The scaling of the transmission values is arbitrary
    because the transmission curve will be normalized by this abstract base class. */
class TabulatedBand : public Band
{
    ITEM_ABSTRACT(TabulatedBand, Band, "a wavelength band (transmission curve) tabulated by the user")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TabulatedBand, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to load the wavelength/transmission pairs and normalizes
        the transmission values as described in the Band class header. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return the wavelengths and the
        corresponding transmission values. The function must guarantee that both arrays have the
        same size. The wavelengths must be listed in increasing or decreasing order. The scaling of
        the transmission values is arbitrary because the transmission curve will be normalized by
        this abstract base class. */
    virtual void getWavelengthsAndTransmissions(Array& lambdav, Array& transv) const = 0;

    //============= Functions required by base class =============

protected:
    /** This function returns the number of elements in the wavelength and transmission data arrays
        held by this subclass. */
    size_t dataSize() const override;

    /** This function returns a pointer to the first wavelength in the corresponding array held by
        this subclass. The number of elements in this array can be obtained through the dataSize()
        function. */
    const double* wavelengthData() const override;

    /** This function returns a pointer to the first transmission value in the corresponding array
        held by this subclass. The number of elements in this array can be obtained through the
        dataSize() function. The transmission values are normalized as described in the Band class
        header. */
    const double* transmissionData() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _lambdav;  // wavelengths
    Array _transv;   // normalized transmission values
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEBAND_HPP
#define FILEBAND_HPP

#include "Band.hpp"

////////////////////////////////////////////////////////////////////

/** A FileBand object represents a wavelength band with a transmission curve that is loaded from an
    input file. The floating point numbers in the first two columns of the text file specify
    respectively the wavelength and the corresponding transmission value. Any additional columns in
    the file are ignored. The transmission outside the range indicated by the first and the last
    wavelength in the list is considered to be zero.

    The wavelengths are by default given in micron (the units can be overridden by column header
    info in the file) and must listed be in increasing order. The transmission values are
    dimensionless. The scaling of these values is arbitrary because the transmission curve will be
    normalized after being loaded. Refer to the description of the Band class for more information.
    */
class FileBand : public Band
{
    ITEM_CONCRETE(FileBand, Band, "a wavelength band (transmission curve) loaded from a text file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(FileBand, "Level2")

    PROPERTY_STRING(filename, "the name of the file defining the transmission curve")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file and normalizes the transmission values as described in
        the Band class header. */
    void setupSelfBefore() override;

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

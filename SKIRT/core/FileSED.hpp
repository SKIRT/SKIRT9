/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILESED_HPP
#define FILESED_HPP

#include "SED.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** A FileSED object represents a spectral energy distribution that is loaded from an input file.
    The floating point numbers in the first two columns of the text file specify respectively the
    wavelength and the corresponding specific luminosity per unit of wavelength. Any additional
    columns in the file are ignored. The wavelengths must be given in micron and must listed be in
    increasing order. The specific luminosity values can be given in arbitrary units because the
    %SED will be normalized after being loaded. The luminosity outside the range indicated by the
    first and the last wavelength in the file is considered to be zero. */
class FileSED : public SED
{
    ITEM_CONCRETE(FileSED, SED, "a spectral energy distribution loaded from a text file")

    PROPERTY_STRING(filename, "the name of the file with the spectral energy distribution")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the input file and precalculates the cumulative distribition for use by
        the other functions in this class. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _inlambdav;   // intrinsic wavelengths (i.e. as read from file)
    Array _inpv;        // intrinsic normalized specific luminosities (i.e. as read from file,
                        //                      but normalized with source range normalization)
    Array _lambdav;     // wavelengths within source range
    Array _pv;          // normalized specific luminosities within source range
    Array _Pv;          // normalized cumulative distribution within source range
};

////////////////////////////////////////////////////////////////////

#endif

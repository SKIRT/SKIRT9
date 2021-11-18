/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABULATEDSED_HPP
#define TABULATEDSED_HPP

#include "Array.hpp"
#include "ContSED.hpp"

////////////////////////////////////////////////////////////////////

/** TabulatedSED is an abstract class for representing spectral energy distributions that are
    tabulated by the user in the form of wavelength/luminosity pairs. The wavelengths must be
    listed in increasing or decreasing order. The luminosity outside the range indicated by the
    first and the last wavelength is considered to be zero.

    The subclass must load the tabulated data, and this abstract class handles everything else. */
class TabulatedSED : public ContSED
{
    ITEM_ABSTRACT(TabulatedSED, ContSED, "a spectral energy distribution tabulated by the user")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TabulatedSED, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function asks the subclass to load the wavelength/luminosity pairs and precalculates
        the cumulative distribution for use by the other functions in this class. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass to return the wavelengths and the
        corresponding luminosities tabulating the %SED. The function must guarantee that both
        arrays have the same size. The wavelengths must be listed in increasing or decreasing
        order. Constant scaling of the luminosities is not important because the %SED will be
        normalized by this abstract class. */
    virtual void getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const = 0;

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For the current class,
        the range is simply obtained from the tabulated wavelengths. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. */
    void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const override;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range. */
    double integratedLuminosity(const Range& wavelengthRange) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution.
        */
    double generateWavelength() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _inlambdav;  // intrinsic wavelengths (i.e. as read from file)
    Array _inpv;       // intrinsic normalized specific luminosities (i.e. as read from file,
                       //                      but normalized with source range normalization)
    Array _lambdav;    // wavelengths within source range
    Array _pv;         // normalized specific luminosities within source range
    Array _Pv;         // normalized cumulative distribution within source range
};

////////////////////////////////////////////////////////////////////

#endif

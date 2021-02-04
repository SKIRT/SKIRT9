/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONTSED_HPP
#define CONTSED_HPP

#include "Array.hpp"
#include "SED.hpp"

////////////////////////////////////////////////////////////////////

/** ContSED is an abstract class for representing regular continuous spectral energy distributions
    that define a specific luminosity at any wavelength in the intrinsic wavelength range. In
    addition to the interface required by the SED base class, ContSED subclasses must implement
    functions to obtain the specific luminosity at a given wavelength or for a sampled range of
    wavelengths. */
class ContSED : public SED
{
    ITEM_ABSTRACT(ContSED, SED, "a continuous spectral energy distribution")
    ITEM_END()

    //============== Functions to be implemented in subclasses ==============

public:
    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength, or zero if the wavelength is
        outside of the %SED's intrinsic wavelength range. */
    virtual double specificLuminosity(double wavelength) const = 0;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. The number of points returned is implementation-dependent and usually
        matches the internal tabular representation of the %SED. The minimum and maximum
        wavelengths in the specified range may or may not be included in the returned result if
        they fall outside of the %SED's intrinsic wavelength range. In that case, the specific
        luminosities outside the returned range should be assumed to be zero. */
    virtual void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

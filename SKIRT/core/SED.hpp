/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SED_HPP
#define SED_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
#include "Range.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** An instance of a SED subclass represents a spectral energy distribution \f$L_\lambda\f$, i.e.
    power per unit of wavelength. This abstract base class just defines an interface that must be
    implemented by each subclass. There are three key operations: drawing a random wavelength from
    the spectral energy distribution, calculating the specific luminosity at a given wavelength,
    and calculating the integrated luminosity over a given wavelength range.

    The spectral energy distribution is automatically normalized to unity over the source
    wavelength range intersected with the intrinsic wavelength range of the distribution.
    Consequently, the random wavelengths returned by the generateWavelength() function will always
    fall inside that range.

    The functions calculating specific or integrated luminosities use just the intrinsic wavelength
    range of the distribution (i.e. they are not limited by the source range). However, for
    consistency, they use the same normalization as the generateWavelength() function. This
    approach makes it possible for a user to configure a luminosity normalization at a wavelength
    (or over a wavelength range) outside of the wavelength range where the source is actually
    emitting. */
class SED : public SimulationItem
{
    ITEM_ABSTRACT(SED, SimulationItem, "a spectral energy distribution")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at the specified wavelength, or zero if the wavelength is
        outside of the distribution's intrinsic spectral range. */
    virtual double specificLuminosity(double wavelength) const = 0;

    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. radiative
        power per unit of wavelength) at a number of wavelength points within the specified
        wavelength range. The number of points returned is implementation-dependent and usually
        matches the internal tabular representation of the %SED. The minimum and maximum
        wavelengths in the specified range are always included in the returned result, even if they
        fall outside of the distribution's intrinsic spectral range (in which case the
        corresponding specific luminosity is set to zero. */
    virtual void specificLuminosityArray(Array& lambdav, Array& pv, const Range& wavelengthRange) const = 0;

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range, or zero if the range is fully outside of the
        distribution's intrinsic spectral range. */
    virtual double integratedLuminosity(const Range& wavelengthRange) const = 0;

    /** This function draws a random wavelength from the normalized spectral energy distribution
        limited to the external wavelength range. */
    virtual double generateWavelength() const = 0;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif

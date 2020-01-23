/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SED_HPP
#define SED_HPP

#include "Array.hpp"
#include "Range.hpp"
#include "SimulationItem.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** An instance of a SED subclass represents a spectral energy distribution \f$L_\lambda\f$, i.e.
    power per unit of wavelength. This abstract base class just defines an interface that must be
    implemented by each subclass. There are three key operations: drawing a random wavelength from
    the spectral energy distribution, calculating the specific luminosity at a given wavelength,
    and calculating the integrated luminosity over a given wavelength range.

    Each SED subclass must ensure that the implemented spectral energy distribution is normalized
    to unity over the normalization wavelength range, which is defined as the intersection of the
    source wavelength range (determined by the SourceSystem) and the intrinsic wavelength range of
    the implemented distribution (determined by the SED subclass). The wavelengthRange()
    convenience function offered by this abstract base class returns this normalization wavelength
    range.

    Specifically, this means that the random wavelengths returned by the generateWavelength()
    function will always fall inside the normalization range. On the other hand, while the
    functions calculating specific and integrated luminosities use the same normalization, they
    operate across the full intrinsic wavelength range of the %SED without being limited by the
    source wavelength range. This makes it possible for a user to configure a luminosity
    normalization at a wavelength (or over a wavelength range) outside of the wavelength range
    where the sources are actually emitting. */
class SED : public SimulationItem
{
    ITEM_ABSTRACT(SED, SimulationItem, "a spectral energy distribution")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //============== Functions to be implemented in subclasses ==============

public:
    /** This function returns the intrinsic wavelength range of the %SED. Outside this range, all
        luminosities are zero. */
    virtual Range intrinsicWavelengthRange() const = 0;

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

    /** This function returns the normalized integrated luminosity \f$L\f$ (i.e. radiative power)
        over the specified wavelength range, or zero if the range is fully outside of the %SED's
        intrinsic wavelength range. */
    virtual double integratedLuminosity(const Range& wavelengthRange) const = 0;

    /** This function returns a random wavelength drawn from the normalized spectral energy
        distribution limited to the normalization wavelength range. */
    virtual double generateWavelength() const = 0;

    //================== Functions implemented here ==================

public:
    /** This function returns the %SED's normalization wavelength range. This range is defined as
        the intersection of the simulation's source wavelength range (obtained from the simulation
        configuration) and the intrinsic wavelength range of the %SED (obtained through the abtract
        intrinsicWavelengthRange() function which must be impemented in each subclass. */
    Range wavelengthRange() const;

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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RANGEWAVELENGTHDISTRIBUTION_HPP
#define RANGEWAVELENGTHDISTRIBUTION_HPP

#include "Range.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** RangeWavelengthDistribution is an intermediate abstract class that serves as a base class for
    wavelength probability distributions that require the user to configure the intrinsic
    wavelength range. Usually this is because the underlying mathematical function (e.g. uniform,
    logarithmic) has no natural wavelength range of its own.

    The range configured by the user for a RangeWavelengthDistribution object is intersected with
    the wavelength range of the associated source (obtained through the SourceWavelengthRangeInterface).
    As a result, the configured minimum and maximum wavelength values can usually be left to their
    default values (defining a very wide wavelength range). */
class RangeWavelengthDistribution : public WavelengthDistribution
{
    ITEM_ABSTRACT(RangeWavelengthDistribution, WavelengthDistribution,
                  "a wavelength probability distribution with a configured wavelength range")

        PROPERTY_DOUBLE(minWavelength, "the shortest wavelength of the wavelength probability distribution")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(minWavelength, "1 pm")

        PROPERTY_DOUBLE(maxWavelength, "the longest wavelength of the wavelength probability distribution")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(maxWavelength, "1 m")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function determines and caches the wavelength range of the wavelength distribution
        (i.e. the intersection of the external wavelength range and the configured intrinsic
        wavelength range) for use by subclasses. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

protected:
    /** This function returns the wavelength range of the wavelength distribution (i.e. the
        intersection of the external wavelength range and the configured intrinsic wavelength
        range) as a service to subclasses. */
    const Range& range() const { return _range; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Range _range;
};

////////////////////////////////////////////////////////////////////

#endif

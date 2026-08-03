/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BLACKBODYSED_HPP
#define BLACKBODYSED_HPP

#include "Array.hpp"
#include "ContSED.hpp"
class PlanckFunction;

////////////////////////////////////////////////////////////////////

/** BlackBodySED is a class that describes black-body spectral energy distributions, i.e. the
    emission spectra of perfect absorbers which are in thermal equilibrium. Such an %SED is
    characterized by the temperature of the object, and its spectrum is the Planck spectrum.

    Whenever a tabular form of the black body %SED is requested, this class uses a spectral
    resolution of \f$R\triangleq\lambda/\Delta\lambda\geq 1000\f$ (see PlanckFunction::cdf()). */
class BlackBodySED : public ContSED
{
    ITEM_CONCRETE(BlackBodySED, ContSED, "a black-body spectral energy distribution")

        PROPERTY_DOUBLE(temperature, "the black body temperature")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "]0 K")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "5000 K")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a Planck function object corresponding to the user-configured
        temperature and sets up the cumulative distribution that will be used to sample random
        wavelengths. */
    void setupSelfBefore() override;

    /** This function destructs the Planck function object. */
    ~BlackBodySED();

    //======================== Other Functions =======================

public:
    /** This function returns the intrinsic wavelength range of the %SED. For the BlackBodySED, the
        intrinsic range is unlimited, so this function returns a range including all representable
        positive floating point values. */
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
    PlanckFunction* _planck{nullptr};
    Array _lambdav;
    Array _pv;
    Array _Pv;
    double _Ltot{0};
};

////////////////////////////////////////////////////////////////////

#endif

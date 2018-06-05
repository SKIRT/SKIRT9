/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SUNSED_HPP
#define SUNSED_HPP

#include "SED.hpp"
#include "Array.hpp"
#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

/** SunSED represents the spectral energy distribution of the Sun. */
class SunSED : public SED
{
    ITEM_CONCRETE(SunSED, SED, "the spectral energy distribution of the Sun")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function opens the stored table resource tabulating the solar %SED and sets up the
        cumulative distribution that will be used to sample random wavelengths. */
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
    StoredTable<1> _table;
    Array _lambdav;
    Array _cdfv;
    double _Ltot{0};
};

////////////////////////////////////////////////////////////////////

#endif

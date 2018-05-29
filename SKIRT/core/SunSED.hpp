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
    /** This function opens the stored table resource tabulating the solar %SED. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. power per
        unit of wavelength) at the specified wavelength, or zero if the wavelength is outside of
        the distribution's spectral range. */
    double specificLuminosity(double wavelength) const override;

    /** This function draws a random wavelength from the normalized spectral energy distribution
        represented by this object. */
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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SUNSED_HPP
#define SUNSED_HPP

#include "ResourceSED.hpp"

////////////////////////////////////////////////////////////////////

/** SunSED represents the spectral energy distribution of the Sun. The %SED is tabulated over a
    wavelength range from 0.0915 \f$\mu\f$m to 160 \f$\mu\f$m. The source of the data is unknown.
    */
class SunSED : public ResourceSED
{
    ITEM_CONCRETE(SunSED, ResourceSED, "the spectral energy distribution of the Sun")
    ITEM_END()

    //======================== Other Functions =======================

    /** This function returns the name of the stored table resource tabulating the %SED. */
    string resourceName() const override;
};

////////////////////////////////////////////////////////////////////

#endif

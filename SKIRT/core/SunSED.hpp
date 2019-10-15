/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SUNSED_HPP
#define SUNSED_HPP

#include "ResourceSED.hpp"

////////////////////////////////////////////////////////////////////

/** SunSED represents the spectral energy distribution of the Sun. The %SED is tabulated over a
    wavelength range from 0.0915 \f$\mu\mathrm{m}\f$ to 160 \f$\mu\mathrm{m}\f$ with the spectral
    resolution shown in the figure below. The source of the data is unknown.

    \image html SunSed.png
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

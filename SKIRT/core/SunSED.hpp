/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SUNSED_HPP
#define SUNSED_HPP

#include "ResourceSED.hpp"

////////////////////////////////////////////////////////////////////

/** SunSED represents the spectral energy distribution of the Sun. */
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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIUNITS_HPP
#define SIUNITS_HPP

#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/** This class provides a system of standard SI units. For example, the units of length and
    distance are both meter. */
class SIUnits : public Units
{
    ITEM_CONCRETE(SIUnits, Units, "SI units")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STELLARUNITS_HPP
#define STELLARUNITS_HPP

#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/** This class provides a system of units appropriate for a stellar environment. For example, the
    unit of length is 1 AU and the unit of distance is 1 parsec. */
class StellarUnits : public Units
{
    ITEM_CONCRETE(StellarUnits, Units, "stellar units (length in AU, distance in pc)")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

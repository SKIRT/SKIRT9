/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EXTRAGALACTICUNITS_HPP
#define EXTRAGALACTICUNITS_HPP

#include "Units.hpp"

////////////////////////////////////////////////////////////////////

/** This class provides a system of units appropriate for an (extra-)galactic environment. For
    example, the unit of length is 1 pc and the unit of distance is 1 Mpc. */
class ExtragalacticUnits : public Units
{
    ITEM_CONCRETE(ExtragalacticUnits, Units, "extragalactic units (length in pc, distance in Mpc)")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

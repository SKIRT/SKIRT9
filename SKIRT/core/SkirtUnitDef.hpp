/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SKIRTUNITDEF_HPP
#define SKIRTUNITDEF_HPP

#include "UnitDef.hpp"

////////////////////////////////////////////////////////////////////

/** The SkirtUnitDef class defines the units and unit systems used by SKIRT for input/output
    purposes (internally, all quantities in SKIRT are represented in SI units). */
class SkirtUnitDef : public UnitDef
{
public:
    /** The default constructor loads the unit and unit system definitions. */
    SkirtUnitDef();
};

////////////////////////////////////////////////////////////////////

#endif

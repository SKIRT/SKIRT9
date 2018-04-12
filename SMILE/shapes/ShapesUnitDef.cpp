/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ShapesUnitDef.hpp"

////////////////////////////////////////////////////////////////////

ShapesUnitDef::ShapesUnitDef()
{
    addUnit("length", "m", 1.);
    addUnit("length", "cm", 1e-2);
    addUnit("length", "mm", 1e-3);
    addUnit("length", "in", 2.54e-2);
    addUnit("length", "inch", 2.54e-2);

    addUnit("area", "m2", 1.);
    addUnit("area", "cm2", 1e-4);
    addUnit("area", "mm2", 1e-6);
    addUnit("area", "in2", 2.54e-2 * 2.54e-2);

    // add exactly one unit system (with a default unit for each quantity)
    // so that we don't need to define a class hierarchy reflecting the unit systems
    addDefaultUnit("SIUnits", "length", "m");
    addDefaultUnit("SIUnits", "area", "m2");
}

////////////////////////////////////////////////////////////////////

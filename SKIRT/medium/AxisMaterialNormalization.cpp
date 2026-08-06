/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AxisMaterialNormalization.hpp"
#include "FatalError.hpp"
#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

double AxisMaterialNormalization::geometryColumnDensity(const Geometry* geom) const
{
    double geometryColumnDensity = 0.;
    switch (_axis)
    {
        case Axis::X: geometryColumnDensity = geom->SigmaX(); break;
        case Axis::Y: geometryColumnDensity = geom->SigmaY(); break;
        case Axis::Z: geometryColumnDensity = geom->SigmaZ(); break;
    }
    if (geometryColumnDensity <= 0.)
        throw FATALERROR("Can't normalize material for geometry with zero column density along selected axis");

    return geometryColumnDensity;
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpticalDepthMaterialNormalization.hpp"
#include "FatalError.hpp"
#include "Geometry.hpp"
#include "MaterialMix.hpp"

//////////////////////////////////////////////////////////////////////

std::pair<double, double> OpticalDepthMaterialNormalization::numberAndMass(const Geometry* geom,
                                                                           const MaterialMix* mix) const
{
    // get the column density of the geometry along the selected axis
    double geometryColumnDensity = 0.;
    switch(_axis)
    {
    case Axis::X: geometryColumnDensity = geom->SigmaX(); break;
    case Axis::Y: geometryColumnDensity = geom->SigmaY(); break;
    case Axis::Z: geometryColumnDensity = geom->SigmaZ(); break;
    }
    if (geometryColumnDensity <= 0.)
        throw FATALERROR("Can't normalize material for geometry with zero column density along selected axis");

    // calculate the requested number and mass column densities from the configured optical depth
    double requestedNumberColumnDensity = _opticalDepth / mix->sectionExt(_wavelength);
    double requestedMassColumnDensity = requestedNumberColumnDensity * mix->mass();

    return std::make_pair(requestedNumberColumnDensity/geometryColumnDensity,
                          requestedMassColumnDensity/geometryColumnDensity);
}

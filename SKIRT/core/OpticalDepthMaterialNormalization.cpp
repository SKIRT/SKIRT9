/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpticalDepthMaterialNormalization.hpp"
#include "Geometry.hpp"
#include "MaterialMix.hpp"

//////////////////////////////////////////////////////////////////////

std::pair<double, double> OpticalDepthMaterialNormalization::numberAndMass(const Geometry* geom,
                                                                           const MaterialMix* mix) const
{
    // get the column density of the geometry along the selected axis
    double geomColumnDensity = geometryColumnDensity(geom);

    // calculate the requested number and mass column densities from the configured optical depth
    double reqNumberColumnDensity = opticalDepth() / mix->sectionExt(wavelength());
    double reqMassColumnDensity = reqNumberColumnDensity * mix->mass();

    return std::make_pair(reqNumberColumnDensity/geomColumnDensity, reqMassColumnDensity/geomColumnDensity);
}

//////////////////////////////////////////////////////////////////////

Range OpticalDepthMaterialNormalization::wavelengthRange() const
{
    return Range(wavelength(), wavelength());
}

//////////////////////////////////////////////////////////////////////

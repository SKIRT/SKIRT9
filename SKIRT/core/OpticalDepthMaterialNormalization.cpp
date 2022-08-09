/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpticalDepthMaterialNormalization.hpp"
#include "FatalError.hpp"
#include "Geometry.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

std::pair<double, double> OpticalDepthMaterialNormalization::numberAndMass(const Geometry* geom,
                                                                           const MaterialMix* mix) const
{
    // get the column density of the geometry along the selected axis
    double geomColumnDensity = geometryColumnDensity(geom);

    // get the extinction cross section of the material at the configured wavelength
    double section = mix->sectionExt(wavelength());
    if (section <= 0.) throw FATALERROR("Can't normalize optical depth for material with zero extinction");

    // calculate the requested number and mass column densities from the configured optical depth
    double reqNumberColumnDensity = opticalDepth() / section;
    double reqMassColumnDensity = reqNumberColumnDensity * mix->mass();

    return std::make_pair(reqNumberColumnDensity / geomColumnDensity, reqMassColumnDensity / geomColumnDensity);
}

////////////////////////////////////////////////////////////////////

Range OpticalDepthMaterialNormalization::wavelengthRange() const
{
    return Range(wavelength(), wavelength());
}

////////////////////////////////////////////////////////////////////

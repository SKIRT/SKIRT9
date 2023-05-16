/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NumberColumnMaterialNormalization.hpp"
#include "Geometry.hpp"
#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

std::pair<double, double> NumberColumnMaterialNormalization::numberAndMass(const Geometry* geom,
                                                                           const MaterialMix* mix) const
{
    // get the column density of the geometry along the selected axis
    double geomColumnDensity = geometryColumnDensity(geom);

    // calculate the requested number and mass column densities from the configured optical depth
    double reqNumberColumnDensity = _numberColumnDensity;
    double reqMassColumnDensity = _numberColumnDensity * mix->mass();

    return std::make_pair(reqNumberColumnDensity / geomColumnDensity, reqMassColumnDensity / geomColumnDensity);
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeanPinteBenchmarkDustMix.hpp"

////////////////////////////////////////////////////////////////////

DustMix::ScatteringMode MeanPinteBenchmarkDustMix::scatteringMode() const
{
    switch (scatteringType())
    {
        case ScatteringType::HenyeyGreenstein: return ScatteringMode::HenyeyGreenstein;
        case ScatteringType::MaterialPhaseFunction: return ScatteringMode::MaterialPhaseFunction;
        case ScatteringType::SphericalPolarization: return ScatteringMode::SphericalPolarization;
    }
    return ScatteringMode::HenyeyGreenstein;  // to satisfy gcc compiler
}

////////////////////////////////////////////////////////////////////

bool MeanPinteBenchmarkDustMix::hasPolarizedScattering() const
{
    return scatteringType() == ScatteringType::SphericalPolarization;
}

//////////////////////////////////////////////////////////////////////

string MeanPinteBenchmarkDustMix::resourceNameForOpticalProps() const
{
    return "MeanPinteBenchmarkOpticalProps";
}

//////////////////////////////////////////////////////////////////////

string MeanPinteBenchmarkDustMix::resourceNameForMuellerMatrix() const
{
    return "MeanPinteBenchmarkMuellerMatrix";
}

//////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeanTrustBenchmarkDustMix.hpp"

////////////////////////////////////////////////////////////////////

MaterialMix::ScatteringMode MeanTrustBenchmarkDustMix::scatteringMode() const
{
    switch (scatteringType())
    {
        case ScatteringType::HenyeyGreenstein: return ScatteringMode::HenyeyGreenstein;
        case ScatteringType::MaterialPhaseFunction: return ScatteringMode::MaterialPhaseFunction;
        case ScatteringType::SphericalPolarization: return ScatteringMode::SphericalPolarization;
    }
    return ScatteringMode::HenyeyGreenstein;  // to satisfy gcc compiler
}

//////////////////////////////////////////////////////////////////////

string MeanTrustBenchmarkDustMix::resourceNameForOpticalProps() const
{
    return "MeanTrustBenchmarkOpticalProps";
}

//////////////////////////////////////////////////////////////////////

string MeanTrustBenchmarkDustMix::resourceNameForMuellerMatrix() const
{
    return "MeanTrustBenchmarkMuellerMatrix";
}

//////////////////////////////////////////////////////////////////////

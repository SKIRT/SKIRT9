/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BlackBodySEDFamily.hpp"
#include "PlanckFunction.hpp"

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> BlackBodySEDFamily::parameterInfo() const
{
    return vector<SnapshotParameter>{
        {"radius", "length", "km"},
        {"temperature", "temperature", "K"},
    };
}

////////////////////////////////////////////////////////////////////

Range BlackBodySEDFamily::intrinsicWavelengthRange() const
{
    return Range(std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::max());
}

////////////////////////////////////////////////////////////////////

double BlackBodySEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double R = parameters[0];
    double T = parameters[1];

    PlanckFunction B(T);
    return 4. * M_PI * M_PI * R * R * B(wavelength);
}

////////////////////////////////////////////////////////////////////

double BlackBodySEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                               const Array& parameters) const
{
    double R = parameters[0];
    double T = parameters[1];

    PlanckFunction B(T);
    return 4. * M_PI * M_PI * R * R * B.cdf(lambdav, pv, Pv, wavelengthRange);
}

////////////////////////////////////////////////////////////////////

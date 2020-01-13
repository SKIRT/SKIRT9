/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileSSPSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

void FileSSPSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, filename(), "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false, false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FileSSPSEDFamily::parameterInfo() const
{
    return vector<SnapshotParameter>{
        {"initial mass", "mass", "Msun"},
        {"metallicity"},
        {"age", "time", "yr"},
    };
}

////////////////////////////////////////////////////////////////////

Range FileSSPSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FileSSPSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double FileSSPSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                             const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////

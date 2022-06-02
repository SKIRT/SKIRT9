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

    if (!_hasIonizationParameter)
        _table3.open(this, filename(), "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false, false);
    else
        _table4.open(this, filename(), "lambda(m),Z(1),t(yr),U(1)", "Llambda(W/m)", false, false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FileSSPSEDFamily::parameterInfo() const
{
    vector<SnapshotParameter> info{SnapshotParameter::initialMass(), SnapshotParameter::metallicity(),
                                   SnapshotParameter::age()};
    if (_hasIonizationParameter) info.push_back(SnapshotParameter::custom("ionization parameter"));
    return info;
}

////////////////////////////////////////////////////////////////////

Range FileSSPSEDFamily::intrinsicWavelengthRange() const
{
    if (!_hasIonizationParameter)
        return _table3.axisRange<0>();
    else
        return _table4.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FileSSPSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    if (!_hasIonizationParameter)
    {
        return M * _table3(wavelength, Z, t);
    }
    else
    {
        double U = parameters[3];
        return M * _table4(wavelength, Z, t, U);
    }
}

////////////////////////////////////////////////////////////////////

double FileSSPSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                             const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    if (!_hasIonizationParameter)
    {
        return M * _table3.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
    }
    else
    {
        double U = parameters[3];
        return M * _table4.cdf(lambdav, pv, Pv, wavelengthRange, Z, t, U);
    }
}

////////////////////////////////////////////////////////////////////

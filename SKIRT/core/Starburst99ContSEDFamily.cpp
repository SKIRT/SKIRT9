/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Starburst99ContSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

Starburst99ContSEDFamily::Starburst99ContSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void Starburst99ContSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "Starburst99ContSEDFamily", "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> Starburst99ContSEDFamily::parameterInfo() const
{
    return vector<SnapshotParameter>{
        {"star formation rate", "massrate", "Msun/yr"},
        {"metallicity"},
        {"age", "time", "yr"},
    };
}

////////////////////////////////////////////////////////////////////

Range Starburst99ContSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double Starburst99ContSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double SFR = parameters[0] / Constants::Msun() * Constants::year();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return SFR * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double Starburst99ContSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                 const Array& parameters) const
{
    double SFR = parameters[0] / Constants::Msun() * Constants::year();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return SFR * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////

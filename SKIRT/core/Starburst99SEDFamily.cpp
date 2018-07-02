/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Starburst99SEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

Starburst99SEDFamily::Starburst99SEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void Starburst99SEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "Starburst99SEDFamily", "lambda(m),Z(1),t(yr)", "Llambda(W/m)");
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> Starburst99SEDFamily::parameterInfo() const
{
    return vector<SnapshotParameter>
    {
        { "initial mass", "mass", "Msun" },
        { "metallicity" },
        { "age", "time", "yr" },
    };
}

////////////////////////////////////////////////////////////////////

double Starburst99SEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double Starburst99SEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv,
                                 const Range& wavelengthRange, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////

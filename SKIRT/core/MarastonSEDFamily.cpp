/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MarastonSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

MarastonSEDFamily::MarastonSEDFamily(SimulationItem* parent, IMF imf)
{
    parent->addChild(this);
    _imf = imf;
    setup();
}

////////////////////////////////////////////////////////////////////

void MarastonSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    string name = "MarastonSEDFamily_";
    name += _imf == IMF::Kroupa ? "Kroupa" : "Salpeter";

    _table.open(this, name, "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> MarastonSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::age()};
}

////////////////////////////////////////////////////////////////////

Range MarastonSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double MarastonSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    // if needed, force the parameter values inside the valid portion of grid
    if ((Z < 0.000894 || Z > 0.0447) && t < 1e9) t = 1e9;

    return M * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double MarastonSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                              const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    // if needed, force the parameter values inside the valid portion of grid
    if ((Z < 0.000894 || Z > 0.0447) && t < 1e9) t = 1e9;

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////

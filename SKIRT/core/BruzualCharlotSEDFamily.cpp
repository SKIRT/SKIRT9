/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BruzualCharlotSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

BruzualCharlotSEDFamily::BruzualCharlotSEDFamily(SimulationItem* parent, IMF imf, Resolution resolution)
{
    parent->addChild(this);
    _imf = imf;
    _resolution = resolution;
    setup();
}

////////////////////////////////////////////////////////////////////

void BruzualCharlotSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    string name = "BruzualCharlotSEDFamily_";
    name += _imf == IMF::Chabrier ? "Chabrier" : "Salpeter";
    name += "_";
    name += _resolution == Resolution::Low ? "lr" : "hr";

    _table.open(this, name, "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> BruzualCharlotSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::age()};
}

////////////////////////////////////////////////////////////////////

Range BruzualCharlotSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double BruzualCharlotSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double BruzualCharlotSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                    const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////

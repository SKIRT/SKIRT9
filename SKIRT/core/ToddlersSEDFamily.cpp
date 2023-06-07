/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

ToddlersSEDFamily::ToddlersSEDFamily(SimulationItem* parent, PAHFraction pahfraction, Resolution resolution)
{
    parent->addChild(this);
    _pahfraction = pahfraction;
    _resolution = resolution;
    setup();
}

////////////////////////////////////////////////////////////////////

void ToddlersSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    string name = "ToddlersSEDFamily_";
    name += _pahfraction == PAHFraction::High ? "high" : "low";
    name += "PAHfrac_";
    name += _resolution == Resolution::Low ? "lr" : "hr";

    _table.open(this, name, "lambda(m),t(Myr),Z(1),SFE(1),n_cl(1/cm3)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> ToddlersSEDFamily::parameterInfo() const
{
    return {
        SnapshotParameter::age(),
        SnapshotParameter::metallicity(),
        SnapshotParameter::custom("Star formation efficiency"),
        SnapshotParameter::custom("Cloud number density", "numbervolumedensity", "1/cm3"),
        SnapshotParameter::custom("Mass", "mass", "Msun"),
    };
}

////////////////////////////////////////////////////////////////////

Range ToddlersSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double ToddlersSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const

{
    double age = parameters[0] / (1e6 * Constants::year());
    double Z = parameters[1];
    double SFE = parameters[2];
    double n_cl = parameters[3] / 1e6;
    double M = parameters[4] / Constants::Msun();

    return M * _table(wavelength, age, Z, SFE, n_cl);
}

////////////////////////////////////////////////////////////////////

double ToddlersSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                              const Array& parameters) const

{
    double age = parameters[0] / (1e6 * Constants::year());
    double Z = parameters[1];
    double SFE = parameters[2];
    double n_cl = parameters[3] / 1e6;
    double M = parameters[4] / Constants::Msun();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, age, Z, SFE, n_cl);
}
////////////////////////////////////////////////////////////////////

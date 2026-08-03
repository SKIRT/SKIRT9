/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

FSPSSEDFamily::FSPSSEDFamily(SimulationItem* parent, IMF imf)
{
    parent->addChild(this);
    _imf = imf;
    setup();
}

////////////////////////////////////////////////////////////////////

void FSPSSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    string name = "FSPSSEDFamily_";
    switch (_imf)
    {
        case IMF::Chabrier: name += "Chabrier"; break;
        case IMF::Kroupa: name += "Kroupa"; break;
        case IMF::Salpeter: name += "Salpeter"; break;
    }

    _table.open(this, name, "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FSPSSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::age()};
}

////////////////////////////////////////////////////////////////////

Range FSPSSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FSPSSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double FSPSSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                          const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////

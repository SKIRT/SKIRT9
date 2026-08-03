/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MappingsSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

MappingsSEDFamily::MappingsSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void MappingsSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "MappingsSEDFamily", "lambda(m),Z(1),logC(1),P(Pa),fPDR(1)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> MappingsSEDFamily::parameterInfo() const
{
    return {
        SnapshotParameter::custom("star formation rate", "massrate", "Msun/yr"),
        SnapshotParameter::metallicity(),
        SnapshotParameter::custom("compactness"),
        SnapshotParameter::custom("pressure", "pressure", "Pa"),
        SnapshotParameter::custom("covering factor"),
    };
}

////////////////////////////////////////////////////////////////////

Range MappingsSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double MappingsSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double SFR = parameters[0] / Constants::Msun() * Constants::year();
    double Z = parameters[1];
    double logC = parameters[2];
    double p = parameters[3];
    double fPDR = parameters[4];

    return SFR * _table(wavelength, Z, logC, p, fPDR);
}

////////////////////////////////////////////////////////////////////

double MappingsSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                              const Array& parameters) const
{
    double SFR = parameters[0] / Constants::Msun() * Constants::year();
    double Z = parameters[1];
    double logC = parameters[2];
    double p = parameters[3];
    double fPDR = parameters[4];

    return SFR * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, logC, p, fPDR);
}

////////////////////////////////////////////////////////////////////

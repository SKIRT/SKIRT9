/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BpassSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

BpassSEDFamily::BpassSEDFamily(SimulationItem* parent, IMF imf, Resolution resolution)
{
    parent->addChild(this);
    _imf = imf;
    _resolution = resolution;
    setup();
}

////////////////////////////////////////////////////////////////////

void BpassSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    string name = "BpassSEDFamily";
    switch (_imf)
    {
        case IMF::Chabrier100: name += "_Chabrier100"; break;
        case IMF::Chabrier300: name += "_Chabrier300"; break;
    }
    if (_resolution == Resolution::Downsampled) name += "_downsampled";

    _table.open(this, name, "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> BpassSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::age()};
}

////////////////////////////////////////////////////////////////////

Range BpassSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double BpassSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double BpassSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                           const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////

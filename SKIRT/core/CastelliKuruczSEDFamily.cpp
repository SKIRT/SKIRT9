/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CastelliKuruczSEDFamily.hpp"

///////////////////////////////////////////////////////////////////

CastelliKuruczSEDFamily::CastelliKuruczSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void CastelliKuruczSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();
    _table.open(this, "CastelliKuruczSEDFamily", "lambda(m),Z(1),Teff(K),g(m/s2)", "Flambda(W/m2/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> CastelliKuruczSEDFamily::parameterInfo() const
{
    return vector<SnapshotParameter>
    {
        { "radius", "length", "km" },
        { "metallicity" },
        { "effective temperature", "temperature", "K" },
        { "surface gravity", "acceleration", "m/s2" },
    };
}

////////////////////////////////////////////////////////////////////

double CastelliKuruczSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double R = parameters[0];
    double Z = parameters[1];
    double T = parameters[2];
    double g = parameters[3];

    // if needed, force the parameter values inside the valid portion of the grid

    return 4.*M_PI * R*R * _table(wavelength, Z, T, g);
}

////////////////////////////////////////////////////////////////////

double CastelliKuruczSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv,
                                    const Range& wavelengthRange, const Array& parameters) const
{
    double R = parameters[0];
    double Z = parameters[1];
    double T = parameters[2];
    double g = parameters[3];

    // if needed, force the parameter values inside the valid portion of the grid

    return 4.*M_PI * R*R * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, T, g);
}

////////////////////////////////////////////////////////////////////

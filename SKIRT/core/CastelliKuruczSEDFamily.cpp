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
    return {SnapshotParameter::custom("radius", "length", "km"), SnapshotParameter::metallicity(),
            SnapshotParameter::temperature(), SnapshotParameter::custom("surface gravity", "acceleration", "m/s2")};
}

////////////////////////////////////////////////////////////////////

Range CastelliKuruczSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // force the gravity value inside the valid portion of the grid depending on the temperature
    void clampParameterValues(double& T, double& g)
    {
        // cutoff values for temperature and gravity (see table in class documentation)
        static Array Tv = {49000, 39000, 31000, 26000, 19000, 11750, 9000, 8250, 7500, 6000};
        static Array gv = pow(10., Array({5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.5}) - 2.);
        static size_t n = Tv.size();

        // clamp the gravity value if needed
        for (size_t i = 0; i != n; ++i)
        {
            if (T > Tv[i] && g < gv[i])
            {
                g = gv[i];
                break;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

double CastelliKuruczSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double R = parameters[0];
    double Z = parameters[1];
    double T = parameters[2];
    double g = parameters[3];

    // if needed, force the parameter values inside the valid portion of the grid
    clampParameterValues(T, g);

    return 4. * M_PI * R * R * _table(wavelength, Z, T, g);
}

////////////////////////////////////////////////////////////////////

double CastelliKuruczSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                    const Array& parameters) const
{
    double R = parameters[0];
    double Z = parameters[1];
    double T = parameters[2];
    double g = parameters[3];

    // if needed, force the parameter values inside the valid portion of the grid
    clampParameterValues(T, g);

    return 4. * M_PI * R * R * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, T, g);
}

////////////////////////////////////////////////////////////////////

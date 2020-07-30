/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileEmissionLineSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

void FileEmissionLineSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, filename(), "lambda(m),index(1)", "Llambda(W/m)", false, false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FileEmissionLineSEDFamily::parameterInfo() const
{
    return vector<SnapshotParameter>{
        {"index"},
    };
}

////////////////////////////////////////////////////////////////////

Range FileEmissionLineSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FileEmissionLineSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double index = parameters[0];
    if (index < 0)
    {
        return 0;
    }
    else
    {
        return _table(wavelength, index);
    }
}

////////////////////////////////////////////////////////////////////

double FileEmissionLineSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                      const Array& parameters) const
{
    double index = parameters[0];
    if (index < 0)
    {
        return 0.;
    }
    else
    {
        return _table.cdf(lambdav, pv, Pv, wavelengthRange, index);
    }
}

////////////////////////////////////////////////////////////////////

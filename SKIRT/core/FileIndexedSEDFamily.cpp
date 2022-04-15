/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileIndexedSEDFamily.hpp"

////////////////////////////////////////////////////////////////////

void FileIndexedSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, filename(), "lambda(m),index(1)", "Llambda(W/m)", false, false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FileIndexedSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::custom("index")};
}

////////////////////////////////////////////////////////////////////

Range FileIndexedSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FileIndexedSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double index = parameters[0];
    // ignore negative indices as per the API
    if (index < 0)
    {
        return 0.;
    }
    else
    {
        return _table(wavelength, index);
    }
}

////////////////////////////////////////////////////////////////////

double FileIndexedSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                 const Array& parameters) const
{
    double index = parameters[0];
    // ignore negative indices as per the API
    if (index < 0)
    {
        return 0.;
    }
    else
    {
        // the sanity of the input index is not verified; if the index does not correspond to a value in the file,
        // cdf() will interpolate from existing index values and this interpolation will most likely be wrong
        return _table.cdf(lambdav, pv, Pv, wavelengthRange, index);
    }
}

////////////////////////////////////////////////////////////////////

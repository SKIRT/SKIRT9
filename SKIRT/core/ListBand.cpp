/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ListBand.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void ListBand::setupSelfBefore()
{
    Band::setupSelfBefore();

    // verify number of configured parameters
    if (_wavelengths.size() != _transmissionValues.size())
        throw FATALERROR("Number of listed luminosities does not match number of listed transmission values");
    if (_transmissionValues.size() < 2) throw FATALERROR("Number of listed transmission values is less than 2");

    // calculate the normalization factor for the transmission values
    size_t size = _transmissionValues.size();
    double norm = 0.;
    for (size_t i = 1; i != size; ++i)
    {
        double dlambda = _wavelengths[i] - _wavelengths[i - 1];
        double trans = 0.5 * (_transmissionValues[i - 1] + _transmissionValues[i]);
        norm += trans * dlambda;
    }

    // copy the transmission values from the configuration and normalize them
    _transv = NR::array(_transmissionValues) / norm;
}

////////////////////////////////////////////////////////////////////

size_t ListBand::dataSize() const
{
    return _transv.size();
}

////////////////////////////////////////////////////////////////////

const double* ListBand::wavelengthData() const
{
    return &_wavelengths[0];
}

////////////////////////////////////////////////////////////////////

const double* ListBand::transmissionData() const
{
    return &_transv[0];
}

////////////////////////////////////////////////////////////////////

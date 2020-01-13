/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MeanListDustMix.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

double MeanListDustMix::getDustProperties(Array& lambdav, Array& kappaextv, Array& albedov, Array& asymmparv) const
{
    // verify number of configured properties
    if (_wavelengths.size() != _extinctionCoefficients.size() || _wavelengths.size() != _albedos.size()
        || _wavelengths.size() != _asymmetryParameters.size())
        throw FATALERROR("Number of listed properties does not match number of listed wavelengths");

    // copy the wavelengths and properties from the configuration
    NR::assign(lambdav, _wavelengths);
    NR::assign(kappaextv, _extinctionCoefficients);
    NR::assign(albedov, _albedos);
    NR::assign(asymmparv, _asymmetryParameters);

    return 1.5e-29;  // in kg/H
}

//////////////////////////////////////////////////////////////////////

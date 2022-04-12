/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SelectDustMixFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> SelectDustMixFamily::parameterInfo() const
{
    return {SnapshotParameter::custom("dustmix index")};
}

////////////////////////////////////////////////////////////////////

const MaterialMix* SelectDustMixFamily::mix(const Array& parameters) const
{
    long numMixes = _dustMixes.size();
    long index = max(0L, min(std::lround(parameters[0]), numMixes - 1));
    return _dustMixes[index];
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MaterialMix.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void MaterialMix::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
    _config = find<Configuration>();
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedScattering() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedAbsorption() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasResonantScattering() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasStochasticDustEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasExtraSpecificState() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasSemiDynamicMediumState() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasLineEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasContinuumEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> MaterialMix::parameterInfo() const
{
    return vector<SnapshotParameter>();
}

////////////////////////////////////////////////////////////////////

void MaterialMix::initializeSpecificState(MaterialState* /*state*/, double /*metallicity*/, double /*temperature*/,
                                          const Array& /*params*/) const
{}

////////////////////////////////////////////////////////////////////

double MaterialMix::asymmpar(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::updateSpecificState(MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* MaterialMix::emissionWavelengthGrid() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::emissivity(const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::emissionSpectrum(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::thetaGrid() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::sectionsAbs(double /*lambda*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::sectionsAbspol(double /*lambda*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::lineEmissionCenters() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::lineEmissionMasses() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::lineEmissionSpectrum(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

double MaterialMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

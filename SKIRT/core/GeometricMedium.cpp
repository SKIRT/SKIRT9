/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricMedium.hpp"
#include "Configuration.hpp"

////////////////////////////////////////////////////////////////////

void GeometricMedium::setupSelfAfter()
{
    Medium::setupSelfAfter();

    // determine normalization
    std::tie(_number, _mass) = normalization()->numberAndMass(geometry(), materialMix());
}

////////////////////////////////////////////////////////////////////

int GeometricMedium::dimension() const
{
    int velocityDimension = hasVelocity() ? velocityDistribution()->dimension() : 1;
    int magneticFieldDimension = hasMagneticField() ? magneticFieldDistribution()->dimension() : 1;
    return max({geometry()->dimension(), velocityDimension, magneticFieldDimension});
}

////////////////////////////////////////////////////////////////////

const MaterialMix* GeometricMedium::mix(Position /*bfr*/) const
{
    return materialMix();
}

////////////////////////////////////////////////////////////////////

const MaterialMix* GeometricMedium::mix() const
{
    return materialMix();
}

////////////////////////////////////////////////////////////////////

bool GeometricMedium::hasVariableMix() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool GeometricMedium::hasVelocity() const
{
    if (velocityDistribution() && velocityMagnitude())
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

Vec GeometricMedium::bulkVelocity(Position bfr) const
{
    return hasVelocity() ? velocityMagnitude() * velocityDistribution()->vector(bfr) : Vec();
}

////////////////////////////////////////////////////////////////////

bool GeometricMedium::hasMagneticField() const
{
    return magneticFieldDistribution() && magneticFieldStrength();
}

////////////////////////////////////////////////////////////////////

Vec GeometricMedium::magneticField(Position bfr) const
{
    return hasMagneticField() ? magneticFieldStrength() * magneticFieldDistribution()->vector(bfr) : Vec();
}

////////////////////////////////////////////////////////////////////

bool GeometricMedium::hasMetallicity() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::metallicity(Position /*bfr*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

bool GeometricMedium::hasTemperature() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::temperature(Position /*bfr*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

bool GeometricMedium::hasParameters() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

void GeometricMedium::parameters(Position /*bfr*/, Array& params) const
{
    params.resize(0);
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::numberDensity(Position bfr) const
{
    return _number * geometry()->density(bfr);
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::number() const
{
    return _number;
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::massDensity(Position bfr) const
{
    return _mass * geometry()->density(bfr);
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::opticalDepthX(double lambda) const
{
    return _number * geometry()->SigmaX() * _materialMix->sectionExt(lambda);
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::opticalDepthY(double lambda) const
{
    return _number * geometry()->SigmaY() * _materialMix->sectionExt(lambda);
}

////////////////////////////////////////////////////////////////////

double GeometricMedium::opticalDepthZ(double lambda) const
{
    return _number * geometry()->SigmaZ() * _materialMix->sectionExt(lambda);
}

////////////////////////////////////////////////////////////////////

Position GeometricMedium::generatePosition() const
{
    return geometry()->generatePosition();
}

////////////////////////////////////////////////////////////////////

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
    int velocityDimension = 1;
    if (!find<Configuration>()->oligochromatic())
    {
        if (velocityZ()) velocityDimension = 2;
        if (velocityX() || velocityY()) velocityDimension = 3;
    }
    return max(geometry()->dimension(), velocityDimension);
}

////////////////////////////////////////////////////////////////////

const MaterialMix*GeometricMedium::mix(Position /*bfr*/) const
{
    return materialMix();
}

////////////////////////////////////////////////////////////////////

bool GeometricMedium::hasVelocity() const
{
    return velocityX() || velocityY() || velocityZ();
}

////////////////////////////////////////////////////////////////////

Vec GeometricMedium::bulkVelocity(Position /*bfr*/) const
{
    return Vec(velocityX(), velocityY(), velocityZ());
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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RotateVectorFieldDecorator.hpp"

//////////////////////////////////////////////////////////////////////

void RotateVectorFieldDecorator::setupSelfBefore()
{
    VectorField::setupSelfBefore();

    // cache frequently used values
    double sinalpha = sin(_eulerAlpha);
    double cosalpha = cos(_eulerAlpha);
    double sinbeta = sin(_eulerBeta);
    double cosbeta = cos(_eulerBeta);
    double singamma = sin(_eulerGamma);
    double cosgamma = cos(_eulerGamma);
    _R11 = cosalpha * cosgamma - sinalpha * cosbeta * singamma;
    _R12 = sinalpha * cosgamma + cosalpha * cosbeta * singamma;
    _R13 = sinbeta * singamma;
    _R21 = -cosalpha * singamma - sinalpha * cosbeta * cosgamma;
    _R22 = -sinalpha * singamma + cosalpha * cosbeta * cosgamma;
    _R23 = sinbeta * cosgamma;
    _R31 = sinalpha * sinbeta;
    _R32 = -cosalpha * sinbeta;
    _R33 = cosbeta;
}

////////////////////////////////////////////////////////////////////

int RotateVectorFieldDecorator::dimension() const
{
    return 3;
}

////////////////////////////////////////////////////////////////////

Vec RotateVectorFieldDecorator::vector(Position bfr) const
{
    // derotate the position to the original position, and then rotate the resulting vector back
    return rotate(_vectorField->vector(derotate(bfr)));
}

////////////////////////////////////////////////////////////////////

Vec RotateVectorFieldDecorator::rotate(Vec bfrorig) const
{
    double xorig = bfrorig.x();
    double yorig = bfrorig.y();
    double zorig = bfrorig.z();
    double x = _R11 * xorig + _R12 * yorig + _R13 * zorig;
    double y = _R21 * xorig + _R22 * yorig + _R23 * zorig;
    double z = _R31 * xorig + _R32 * yorig + _R33 * zorig;
    return Vec(x, y, z);
}

////////////////////////////////////////////////////////////////////

Position RotateVectorFieldDecorator::derotate(Position bfr) const
{
    double x = bfr.x();
    double y = bfr.y();
    double z = bfr.z();
    double xorig = _R11 * x + _R21 * y + _R31 * z;
    double yorig = _R12 * x + _R22 * y + _R32 * z;
    double zorig = _R13 * x + _R23 * y + _R33 * z;
    return Position(xorig, yorig, zorig);
}

////////////////////////////////////////////////////////////////////

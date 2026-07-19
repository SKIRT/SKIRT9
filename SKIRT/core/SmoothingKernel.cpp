/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SmoothingKernel.hpp"
#include "Box.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void SmoothingKernel::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // cache random generator
    _random = find<Random>();

    // open the cumulative kernel table used in massInBox()
    _cumkernel.open(this, type(), "X(1),Y(1),Z(1)", "Phi(1)", false);
}

//////////////////////////////////////////////////////////////////////

double SmoothingKernel::massInBox(const Box& box) const
{
    // Phi(X,Y,Z) for positive X,Y,Z: interpolating the table
    auto Phi = [this](double X, double Y, double Z) {
        if (X >= 1.0 && Y >= 1.0 && Z >= 1.0) return 0.125;
        X = std::min(X, 1.0);
        Y = std::min(Y, 1.0);
        Z = std::min(Z, 1.0);
        return _cumkernel(X, Y, Z);
    };

    // Signed Phi*(X,Y,Z) for arbitrary X,Y,Z
    auto PhiSigned = [Phi](double X, double Y, double Z) {
        double s = 1.0;
        if (X < 0.0) s = -s;
        if (Y < 0.0) s = -s;
        if (Z < 0.0) s = -s;
        return s * Phi(std::abs(X), std::abs(Y), std::abs(Z));
    };

    // Performing the integration with 8 values from the cumulative kernel
    double xm, ym, zm, xp, yp, zp;
    box.extent(xm, ym, zm, xp, yp, zp);
    return PhiSigned(xp, yp, zp) - PhiSigned(xp, yp, zm) - PhiSigned(xp, ym, zp) + PhiSigned(xp, ym, zm)
           - PhiSigned(xm, yp, zp) + PhiSigned(xm, yp, zm) + PhiSigned(xm, ym, zp) - PhiSigned(xm, ym, zm);
}

//////////////////////////////////////////////////////////////////////

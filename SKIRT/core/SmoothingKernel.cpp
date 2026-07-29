/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SmoothingKernel.hpp"
#include "Box.hpp"
#include "FilePaths.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void SmoothingKernel::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // cache random generator
    _random = find<Random>();

    // open the cumulative kernel table used in massInBox()
    if (hasMassInBox()) _cumkernel.open(this, type(), "X(1),Y(1),Z(1)", "Phi(1)", false);
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

        // we could simply do: return _cumkernel(X, Y, Z);
        // however, this would cost a binary search along each axis, so we perform the interpolation here,
        // assuming that the three axes have range [0..1] with the same number of equidistant points
        size_t num = _cumkernel.axisSize<0>() - 1;
        X *= num;
        Y *= num;
        Z *= num;
        size_t i0 = std::min(num - 1, static_cast<size_t>(X));
        size_t j0 = std::min(num - 1, static_cast<size_t>(Y));
        size_t k0 = std::min(num - 1, static_cast<size_t>(Z));
        size_t i1 = i0 + 1;
        size_t j1 = j0 + 1;
        size_t k1 = k0 + 1;
        double tx = X - i0;
        double ty = Y - j0;
        double tz = Z - k0;

        double v000 = _cumkernel.valueAtIndices(i0, j0, k0);
        double v100 = _cumkernel.valueAtIndices(i1, j0, k0);
        double v010 = _cumkernel.valueAtIndices(i0, j1, k0);
        double v110 = _cumkernel.valueAtIndices(i1, j1, k0);
        double v001 = _cumkernel.valueAtIndices(i0, j0, k1);
        double v101 = _cumkernel.valueAtIndices(i1, j0, k1);
        double v011 = _cumkernel.valueAtIndices(i0, j1, k1);
        double v111 = _cumkernel.valueAtIndices(i1, j1, k1);
        double v00 = (1.0 - tx) * v000 + tx * v100;
        double v10 = (1.0 - tx) * v010 + tx * v110;
        double v01 = (1.0 - tx) * v001 + tx * v101;
        double v11 = (1.0 - tx) * v011 + tx * v111;
        double v0 = (1.0 - ty) * v00 + ty * v10;
        double v1 = (1.0 - ty) * v01 + ty * v11;
        return (1.0 - tz) * v0 + tz * v1;
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

bool SmoothingKernel::hasMassInBox() const
{
    return FilePaths::hasResource(type() + ".stab");
}

//////////////////////////////////////////////////////////////////////

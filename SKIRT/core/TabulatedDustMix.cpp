/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedDustMix.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

double TabulatedDustMix::getOpticalProperties(const Array& lambdav, const Array& /*thetav*/, Array& sigmaabsv,
                                              Array& sigmascav, Array& asymmparv, Table<2>& /*S11vv*/,
                                              Table<2>& /*S12vv*/, Table<2>& /*S33vv*/, Table<2>& /*S34vv*/,
                                              ArrayTable<2>& /*sigmaabsvv*/, ArrayTable<2>& /*sigmaabspolvv*/)
{
    // obtain the basic properties from the subclass on a wavelength grid chosen by the subclass
    Array inlambdav, inkappaextv, inalbedov, inasymmparv;
    double mu = getDustProperties(inlambdav, inkappaextv, inalbedov, inasymmparv);

    // verify the number of points
    if (inlambdav.size() < 1) throw FATALERROR("Dust properties must be tabulated for at least one wavelength");

    // reverse the arrays if needed to get the wavelengths in increasing order
    if (inlambdav.size() > 1 && inlambdav[0] > inlambdav[inlambdav.size() - 1])
    {
        std::reverse(begin(inlambdav), end(inlambdav));
        std::reverse(begin(inkappaextv), end(inkappaextv));
        std::reverse(begin(inalbedov), end(inalbedov));
        std::reverse(begin(inasymmparv), end(inasymmparv));
    }

    // convert to properties required by caller
    Array insigmaabsv = mu * inkappaextv * (1. - inalbedov);
    Array insigmascav = mu * inkappaextv * inalbedov;

    // resample to wavelength grid required by caller
    sigmaabsv = NR::clampedResample<NR::interpolateLogLog>(lambdav, inlambdav, insigmaabsv);
    sigmascav = NR::clampedResample<NR::interpolateLogLog>(lambdav, inlambdav, insigmascav);
    asymmparv = NR::clampedResample<NR::interpolateLogLin>(lambdav, inlambdav, inasymmparv);

    return mu;
}

//////////////////////////////////////////////////////////////////////

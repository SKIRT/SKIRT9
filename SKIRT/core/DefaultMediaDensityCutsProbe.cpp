/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultMediaDensityCutsProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "PlanarMediaDensityCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

void DefaultMediaDensityCutsProbe::probeSetup()
{
    // The size in pixels for default cuts (in each spatial direction)
    const int Np = 1024;

    if (find<Configuration>()->hasMedium())
    {
        int dimension = find<MediumSystem>()->dimension();

        // output planes depending on the dimension of the medium system
        PlanarMediaDensityCutsProbe::writeMediaDensityCuts(this, 1,1,0, 0.,0.,0., Np,Np,Np);                     // xy
        if (dimension >= 2) PlanarMediaDensityCutsProbe::writeMediaDensityCuts(this, 1,0,1, 0.,0.,0., Np,Np,Np); // xz
        if (dimension == 3) PlanarMediaDensityCutsProbe::writeMediaDensityCuts(this, 0,1,1, 0.,0.,0., Np,Np,Np); // yz
    }
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultCustomStateCutsProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "PlanarCustomStateCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

void DefaultCustomStateCutsProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // the size in pixels (in each spatial direction) for the default cuts
        const int Np = 1024;

        // the dimension of the medium system
        int dimension = find<MediumSystem>()->dimension();

        // output cuts depending on the dimension of the medium system
        PlanarCustomStateCutsProbe::writeCustomStateCuts(this, 1, 1, 0, 0., 0., 0., Np, Np, Np);  // xy
        if (dimension >= 2)
            PlanarCustomStateCutsProbe::writeCustomStateCuts(this, 1, 0, 1, 0., 0., 0., Np, Np, Np);  // xz
        if (dimension == 3)
            PlanarCustomStateCutsProbe::writeCustomStateCuts(this, 0, 1, 1, 0., 0., 0., Np, Np, Np);  // yz
    }
}

////////////////////////////////////////////////////////////////////

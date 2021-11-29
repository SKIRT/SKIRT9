/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultElectronTemperatureCutsProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "PlanarElectronTemperatureCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

void DefaultElectronTemperatureCutsProbe::probeSetup()
{
    if (find<Configuration>()->hasMedium() && find<MediumSystem>()->hasElectrons())
    {
        // the size in pixels (in each spatial direction) for the default cuts
        const int Np = 1024;

        // the dimension of the medium system
        int dimension = find<MediumSystem>()->dimension();

        // output cuts depending on the dimension of the medium system
        PlanarElectronTemperatureCutsProbe::writeElectronTemperatureCut(this, 1, 1, 0, 0., 0., 0., Np, Np, Np);
        if (dimension >= 2)
            PlanarElectronTemperatureCutsProbe::writeElectronTemperatureCut(this, 1, 0, 1, 0., 0., 0., Np, Np, Np);
        if (dimension == 3)
            PlanarElectronTemperatureCutsProbe::writeElectronTemperatureCut(this, 0, 1, 1, 0., 0., 0., Np, Np, Np);
    }
}

////////////////////////////////////////////////////////////////////

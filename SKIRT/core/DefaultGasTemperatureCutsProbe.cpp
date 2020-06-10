/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultGasTemperatureCutsProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "PlanarGasTemperatureCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

void DefaultGasTemperatureCutsProbe::probeSetup()
{
    if (find<Configuration>()->hasLymanAlpha())
    {
        // the size in pixels (in each spatial direction) for the default cuts
        const int Np = 1024;

        // the dimension of the medium system
        int dimension = find<MediumSystem>()->dimension();

        // output cuts depending on the dimension of the medium system
        PlanarGasTemperatureCutsProbe::writeGasTemperatureCut(this, 1, 1, 0, 0., 0., 0., Np, Np, Np);
        if (dimension >= 2)
            PlanarGasTemperatureCutsProbe::writeGasTemperatureCut(this, 1, 0, 1, 0., 0., 0., Np, Np, Np);
        if (dimension == 3)
            PlanarGasTemperatureCutsProbe::writeGasTemperatureCut(this, 0, 1, 1, 0., 0., 0., Np, Np, Np);
    }
}

////////////////////////////////////////////////////////////////////

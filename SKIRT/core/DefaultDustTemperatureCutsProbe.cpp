/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultDustTemperatureCutsProbe.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "PlanarDustTemperatureCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

void DefaultDustTemperatureCutsProbe::probeRun()
{
    if (find<Configuration>()->hasPanRadiationField() && find<MediumSystem>()->hasDust())
    {
        // the size in pixels (in each spatial direction) for the default cuts
        const int Np = 1024;

        // the dimension of the medium system
        int dimension = find<MediumSystem>()->dimension();

        // output cuts depending on the dimension of the medium system
        PlanarDustTemperatureCutsProbe::writeDustTemperatureCut(this, 1,1,0, 0.,0.,0., Np,Np,Np);
        if (dimension >= 2) PlanarDustTemperatureCutsProbe::writeDustTemperatureCut(this, 1,0,1, 0.,0.,0., Np,Np,Np);
        if (dimension == 3) PlanarDustTemperatureCutsProbe::writeDustTemperatureCut(this, 0,1,1, 0.,0.,0., Np,Np,Np);
    }
}

////////////////////////////////////////////////////////////////////

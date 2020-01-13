/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultRadiationFieldCutsProbe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"
#include "PlanarRadiationFieldCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

void DefaultRadiationFieldCutsProbe::probeRun()
{
    if (find<Configuration>()->hasRadiationField())
    {
        // the size in pixels (in each spatial direction) for the default cuts
        const int Np = 1024;

        // the dimension of the medium system
        int dimension = find<MediumSystem>()->dimension();

        // output cuts depending on the dimension of the medium system
        PlanarRadiationFieldCutsProbe::writeRadiationFieldCut(this, 1, 1, 0, 0., 0., 0., Np, Np, Np);
        if (dimension >= 2)
            PlanarRadiationFieldCutsProbe::writeRadiationFieldCut(this, 1, 0, 1, 0., 0., 0., Np, Np, Np);
        if (dimension == 3)
            PlanarRadiationFieldCutsProbe::writeRadiationFieldCut(this, 0, 1, 1, 0., 0., 0., Np, Np, Np);

        // if requested, also output the wavelength grid
        if (writeWavelengthGrid())
        {
            InstrumentWavelengthGridProbe::writeWavelengthGrid(this, find<Configuration>()->radiationFieldWLG(),
                                                               itemName() + "_wavelengths",
                                                               "wavelengths for mean intensity");
        }
    }
}

////////////////////////////////////////////////////////////////////

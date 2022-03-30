/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DefaultPlanarCutsForm.hpp"
#include "MediumSystem.hpp"
#include "PlanarCutsForm.hpp"

////////////////////////////////////////////////////////////////////

void DefaultPlanarCutsForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    // the size in pixels (in each spatial direction) for the default cuts
    constexpr int Np = 1024;

    // the dimension of the medium system
    int dimension = find<MediumSystem>()->dimension();

    // output cuts depending on the dimension of the medium system
    PlanarCutsForm::writePlanarCut(bridge, 1, 1, 0, 0., 0., 0., 0., 0., 0., Np, Np, Np);                      // xy
    if (dimension >= 2) PlanarCutsForm::writePlanarCut(bridge, 1, 0, 1, 0., 0., 0., 0., 0., 0., Np, Np, Np);  // xz
    if (dimension == 3) PlanarCutsForm::writePlanarCut(bridge, 0, 1, 1, 0., 0., 0., 0., 0., 0., Np, Np, Np);  // yz
}

////////////////////////////////////////////////////////////////////

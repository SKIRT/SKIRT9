/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClearDensityRecipe.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

void ClearDensityRecipe::beginUpdate(int /*numCells*/)
{
    _dlambdav = find<Configuration>()->radiationFieldWLG()->dlambdav();
}

////////////////////////////////////////////////////////////////////

UpdateStatus ClearDensityRecipe::update(MaterialState* state, const Array& Jv)
{
    UpdateStatus status;
    if (state->numberDensity() > 0.)
    {
        // the local radiation field in the Milky Way (Mathis et al. 1983) integrated over all wavelengths
        const double JtotMW = 1.7623e-06;
        double U = (Jv * _dlambdav).sum() / JtotMW;
        if (U > fieldStrengthThreshold())
        {
            state->setNumberDensity(0.);
            status.updateNotConverged();
        }
    }
    return status;
}

////////////////////////////////////////////////////////////////////

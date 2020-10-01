/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClearDensityRecipe.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "MediumSystem.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

void ClearDensityRecipe::beginUpdate(int /*numCells*/)
{
    _dlambdav = find<Configuration>()->radiationFieldWLG()->dlambdav();
}

////////////////////////////////////////////////////////////////////

bool ClearDensityRecipe::update(MaterialState* state, const Array& Jv)
{
    if (state->numberDensity() > 0.)
    {
        // the local radiation field in the Milky Way (Mathis et al. 1983) integrated over all wavelengths
        const double JtotMW = 1.7623e-06;
        double U = (Jv * _dlambdav).sum() / JtotMW;
        if (U > fieldStrengthThreshold())
        {
            state->setNumberDensity(0.);
            return true;
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////

bool ClearDensityRecipe::endUpdate(int numCells, int numUpdated)
{
    if (numUpdated == 0)
    {
        find<Log>()->info("ClearDensityRecipe has converged: no cells have been cleared");
        return true;
    }
    else
    {
        double percentage = 100. * numUpdated / numCells;
        find<Log>()->info("ClearDensityRecipe has not yet converged: " + std::to_string(numUpdated) + " out of "
                          + std::to_string(numCells) + " cells (" + StringUtils::toString(percentage, 'f', 2)
                          + " %) have been cleared");
        return false;
    }
}

////////////////////////////////////////////////////////////////////

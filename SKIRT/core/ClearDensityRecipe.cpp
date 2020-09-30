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

void ClearDensityRecipe::beginUpdate()
{
    _dlambdav = find<Configuration>()->radiationFieldWLG()->dlambdav();
    _numCleared = 0;
}

////////////////////////////////////////////////////////////////////

void ClearDensityRecipe::update(MaterialState* state, const Array& Jv)
{
    if (state->numberDensity() > 0.)
    {
        // the local radiation field in the Milky Way (Mathis et al. 1983) integrated over all wavelengths
        const double JtotMW = 1.7623e-06;
        double U = (Jv * _dlambdav).sum() / JtotMW;
        if (U > fieldStrengthThreshold())
        {
            state->setNumberDensity(0.);
            ++_numCleared;
        }
    }
}

////////////////////////////////////////////////////////////////////

bool ClearDensityRecipe::endUpdate()
{
    if (_numCleared == 0)
    {
        find<Log>()->info("ClearDensityRecipe has converged: no cells have been cleared");
        return true;
    }
    else
    {
        auto ms = find<MediumSystem>();
        int numClearedCells = _numCleared / ms->numMedia();
        double percentage = 100. * numClearedCells / ms->numCells();
        find<Log>()->info("ClearDensityRecipe has not yet converged: " + std::to_string(numClearedCells) + " out of "
                   + std::to_string(ms->numCells()) + " cells (" + StringUtils::toString(percentage, 'f', 2)
                   + " %) have been cleared");
        return false;
    }
}

////////////////////////////////////////////////////////////////////

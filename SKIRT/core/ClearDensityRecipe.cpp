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

void ClearDensityRecipe::setupSelfBefore()
{
    DynamicStateRecipe::setupSelfBefore();

    _log = find<Log>();
    _dlambdav = find<Configuration>()->radiationFieldWLG()->dlambdav();

    auto ms = find<MediumSystem>();
    _numCells = ms->numCells();
    _numMedia = ms->numMedia();
}

//////////////////////////////////////////////////////////////////////

void ClearDensityRecipe::beginUpdate()
{
    _numCleared = 0;
}

//////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////

bool ClearDensityRecipe::endUpdate()
{
    if (_numCleared == 0)
    {
        _log->info("ClearDensityRecipe has converged: no cells have been cleared");
        return true;
    }
    else
    {
        int numClearedCells = _numCleared / _numMedia;
        double percentage = 100. * numClearedCells / _numCells;
        _log->info("ClearDensityRecipe has not yet converged: " + std::to_string(numClearedCells) + " out of "
                   + std::to_string(_numCells) + " cells (" + StringUtils::toString(percentage, 'f', 3)
                   + " %) have been cleared)");
        return false;
    }
}

//////////////////////////////////////////////////////////////////////

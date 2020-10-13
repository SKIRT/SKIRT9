/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClearDensityRecipe.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool DynamicStateRecipe::endUpdate(int numCells, int numUpdated)
{
    if (numUpdated == 0)
    {
        find<Log>()->info(type() + " has converged: no cells have been updated");
        return true;
    }
    else
    {
        double percentage = 100. * numUpdated / numCells;
        find<Log>()->info(type() + " has not yet converged: " + std::to_string(numUpdated) + " out of "
                          + std::to_string(numCells) + " cells (" + StringUtils::toString(percentage, 'f', 2)
                          + " %) have been updated");
        return false;
    }
}

////////////////////////////////////////////////////////////////////

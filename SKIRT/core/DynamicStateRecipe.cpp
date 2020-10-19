/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClearDensityRecipe.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool DynamicStateRecipe::endUpdate(int numCells, int numUpdated, int numNotConverged)
{
    // log converged or not
    if (numNotConverged <= maxNotConvergedCells())
        find<Log>()->info(type() + " has converged");
    else
        find<Log>()->info(type() + " has not yet converged");

    // log extra info
    find<Log>()->info("  Updated cells: " + std::to_string(numUpdated) + " out of " + std::to_string(numCells)
                      + " (" + StringUtils::toString(100. * numUpdated / numCells, 'f', 2) + " %)");
    find<Log>()->info("  Not converged: " + std::to_string(numNotConverged) + " out of " + std::to_string(numCells)
                      + " (" + StringUtils::toString(100. * numNotConverged / numCells, 'f', 2) + " %)");

    return numNotConverged <= maxNotConvergedCells();
}

////////////////////////////////////////////////////////////////////

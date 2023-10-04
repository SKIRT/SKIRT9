/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NestedDensityTreePolicy.hpp"
#include "Log.hpp"
#include "SpatialGrid.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

void NestedDensityTreePolicy::setupSelfBefore()
{
    DensityTreePolicy::setupSelfBefore();

    _inner = Box(_innerMinX, _innerMinY, _innerMinZ, _innerMaxX, _innerMaxY, _innerMaxZ);

    auto grid = find<SpatialGrid>(false);
    if (!grid->boundingBox().contains(_inner))
        find<Log>()->warning("The nested density tree policy region is not fully inside the spatial grid domain");
}

////////////////////////////////////////////////////////////////////

bool NestedDensityTreePolicy::needsSubdivide(TreeNode* node)
{
    if (_inner.intersects(node->extent()))
        return _innerPolicy->needsSubdivide(node);
    else
        return DensityTreePolicy::needsSubdivide(node);
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NestedDensityTreePolicy.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "SpatialGrid.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

void NestedDensityTreePolicy::setupSelfBefore()
{
    DensityTreePolicy::setupSelfBefore();

    _inner = Box(_minXInner, _minYInner, _minZInner, _maxXInner, _maxYInner, _maxZInner);

    auto ms = find<MediumSystem>(false);
    if (!ms->grid()->boundingBox().contains(_inner))
        find<Log>()->warning("The nested density tree policy region is not fully inside the spatial grid domain");
}

////////////////////////////////////////////////////////////////////

bool NestedDensityTreePolicy::needsSubdivide(TreeNode* node)
{
    if (_inner.intersects(node->extent()))
        return _nestedTree->needsSubdivide(node);
    else
        return DensityTreePolicy::needsSubdivide(node);
}

////////////////////////////////////////////////////////////////////

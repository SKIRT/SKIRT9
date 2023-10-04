/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NestedDensityTreePolicy.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "SpatialGrid.hpp"
#include "TreeNode.hpp"

////////////////////////////////////////////////////////////////////

void NestedDensityTreePolicy::setupSelfBefore()
{
    DensityTreePolicy::setupSelfBefore();

    // verify that the inner box is not empty
    if (_innerMaxX <= _innerMinX) throw FATALERROR("The extent of the inner box should be positive in the X direction");
    if (_innerMaxY <= _innerMinY) throw FATALERROR("The extent of the inner box should be positive in the Y direction");
    if (_innerMaxZ <= _innerMinZ) throw FATALERROR("The extent of the inner box should be positive in the Z direction");

    // copy the coordinates to a Box instance for ease of use
    _inner = Box(_innerMinX, _innerMinY, _innerMinZ, _innerMaxX, _innerMaxY, _innerMaxZ);

    // verify that the inner box is inside the spatial grid
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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NestedDensityTreePolicy.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "TreeNode.hpp"
#include "TreeSpatialGrid.hpp"

void NestedDensityTreePolicy::setupSelfBefore()
{
    DensityTreePolicy::setupSelfBefore();

    _inner = Box(_minXInner, _minYInner, _minZInner, _maxXInner, _maxYInner, _maxZInner);

    auto ms = find<MediumSystem>(false);
    auto grid = ms->find<TreeSpatialGrid>(false);
    Box _outer = grid->boundingBox();

    if (!_outer.contains(_inner)) throw FATALERROR("The nested density tree is not contained in the outer region");

    _outerMinLevel = _minLevel;
    _outerMaxLevel = _maxLevel;
}

void NestedDensityTreePolicy::setupSelfAfter()
{
    DensityTreePolicy::setupSelfAfter();

    _minLevel = min(_minLevel, _nestedTree->_minLevel);
    _maxLevel = max(_maxLevel, _nestedTree->_maxLevel);
}

bool NestedDensityTreePolicy::needsSubdivide(TreeNode* node, int level)
{
    if (_inner.intersects(node->extent()))
    {
        if (level < _nestedTree->_minLevel) return true;
        if (level >= _nestedTree->_maxLevel) return false;

        return _nestedTree->needsSubdivide(node, level);
    }
    else
    {
        if (level < _outerMinLevel) return true;
        if (level >= _outerMaxLevel) return false;

        return DensityTreePolicy::needsSubdivide(node, level);
    }
}
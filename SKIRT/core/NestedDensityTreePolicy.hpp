/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NESTEDDENSITYTREEPOLICY_HPP
#define NESTEDDENSITYTREEPOLICY_HPP

#include "Box.hpp"
#include "DensityTreePolicy.hpp"

/** NestedDensityTreePolicy is a tree policy that allows for nesting DensityTreePolicies, enabling
    users to define a region containing another (nested) density tree.

    This class inherits from DensityTreePolicy and uses all the inherited properties to classify the
    outer region. Additionally, it introduces an additional property called nestedTree, which is
    also a DensityTreePolicy. The nestedTree property decides the properties of the inner region.

    If a TreeNode intersects with the inner 'nested' region the node will be subdivided based on the
    properties of the nestedTree. This means that even TreeNodes almost fully outside the inner
    region can be subdivided if they have a non-zero intersection with the inner region. */
class NestedDensityTreePolicy : public DensityTreePolicy
{
    ITEM_CONCRETE(NestedDensityTreePolicy, DensityTreePolicy,
                  "a tree grid construction policy using a nested density tree policy")
        ATTRIBUTE_TYPE_DISPLAYED_IF(NestedDensityTreePolicy, "Level2")

        PROPERTY_DOUBLE(minXInner, "the start point of the inner box in the X direction")
        ATTRIBUTE_QUANTITY(minXInner, "length")

        PROPERTY_DOUBLE(maxXInner, "the end point of the inner box in the X direction")
        ATTRIBUTE_QUANTITY(maxXInner, "length")

        PROPERTY_DOUBLE(minYInner, "the start point of the inner box in the Y direction")
        ATTRIBUTE_QUANTITY(minYInner, "length")

        PROPERTY_DOUBLE(maxYInner, "the end point of the inner box in the Y direction")
        ATTRIBUTE_QUANTITY(maxYInner, "length")

        PROPERTY_DOUBLE(minZInner, "the start point of the inner box in the Z direction")
        ATTRIBUTE_QUANTITY(minZInner, "length")

        PROPERTY_DOUBLE(maxZInner, "the end point of the inner box in the Z direction")
        ATTRIBUTE_QUANTITY(maxZInner, "length")

        PROPERTY_ITEM(nestedTree, DensityTreePolicy, "the nested density tree inside the inner region")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function checks whether the inner region is inside the outer region. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function determines whether the given node needs to be subdivided. Depending on
        whether the node intersects with the inner region or not, the request is passed to the
        needsSubdivide() function of the nested policy or to the needsSubdivide() function of our
        base class. */
    bool needsSubdivide(TreeNode* node) override;

    //======================== Data Members ========================

private:
    // the inner region box
    Box _inner;
};

#endif

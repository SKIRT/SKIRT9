/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NESTEDDENSITYTREEPOLICY_HPP
#define NESTEDDENSITYTREEPOLICY_HPP

#include "Box.hpp"
#include "DensityTreePolicy.hpp"

/** NestedDensityTreePolicy is a DensityTreePolicy policy that allows defining separate subdivision
    criteria in a given subregion of the spatial grid. This can be used, for example, to specify a
    higher resolution in a given region of interest.

    The class inherits from DensityTreePolicy and uses all the inherited properties to classify the
    outer region. Additional \em innerMin/innerMax coordinate properties define the bounding box of
    the inner region. Finally, the additional property \em innerPolicy, which is also expected to
    be a DensityTreePolicy instance, specifies the subdivision criteria for the inner region.

    If a TreeNode intersects with the inner region the node will be subdivided based on the
    criteria defined by the \em innerPolicy. This means that even TreeNodes that are almost fully
    outside the inner region can be subdivided by those criteria if they have a non-zero
    intersection with the inner region.

    It is not meaningful to specify an inner region that extends outside of the spatial grid
    bounding box, because none of the tree nodes will intersect those outside areas. Therefore, a
    warning is issued if the inner region is not fully inside the spatial grid domain.

    <b>Recursive nesting</b>

    By default, the inner policy is a regular DensityTreePolicy instance. However, because
    NestedDensityTreePolicy inherits DensityTreePolicy, it is also possible to again select a
    NestedDensityTreePolicy instance as the inner policy, leading to recursive nesting. While this
    is not a recommended use case, it would allow specifying recursively increasing resolution in
    nested, successively smaller regions of the spatial domain. */
class NestedDensityTreePolicy : public DensityTreePolicy
{
    ITEM_CONCRETE(NestedDensityTreePolicy, DensityTreePolicy,
                  "a tree grid construction policy using a nested density tree policy")
        ATTRIBUTE_TYPE_DISPLAYED_IF(NestedDensityTreePolicy, "Level2")

        PROPERTY_DOUBLE(innerMinX, "the start point of the inner box in the X direction")
        ATTRIBUTE_QUANTITY(innerMinX, "length")

        PROPERTY_DOUBLE(innerMaxX, "the end point of the inner box in the X direction")
        ATTRIBUTE_QUANTITY(innerMaxX, "length")

        PROPERTY_DOUBLE(innerMinY, "the start point of the inner box in the Y direction")
        ATTRIBUTE_QUANTITY(innerMinY, "length")

        PROPERTY_DOUBLE(innerMaxY, "the end point of the inner box in the Y direction")
        ATTRIBUTE_QUANTITY(innerMaxY, "length")

        PROPERTY_DOUBLE(innerMinZ, "the start point of the inner box in the Z direction")
        ATTRIBUTE_QUANTITY(innerMinZ, "length")

        PROPERTY_DOUBLE(innerMaxZ, "the end point of the inner box in the Z direction")
        ATTRIBUTE_QUANTITY(innerMaxZ, "length")

        PROPERTY_ITEM(innerPolicy, DensityTreePolicy, "the density tree policy for the inner region")
        ATTRIBUTE_DEFAULT_VALUE(innerPolicy, "DensityTreePolicy")

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

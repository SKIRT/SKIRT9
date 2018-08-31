/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DENSITYTREEPOLICY_HPP
#define DENSITYTREEPOLICY_HPP

#include "TreePolicy.hpp"

//////////////////////////////////////////////////////////////////////

/** DensityTreePolicy represents the configurable options and the corresponding implementation
    mechanisms for constructing spatial tree grids based on the density distribution of the media
    in the medium system.

    TO DO. */
class DensityTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(DensityTreePolicy, TreePolicy,
                  "a tree grid construction policy using the medium density distribution")

    // TO DO.

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function constructs the hierarchical tree and all (interconnected) nodes forming the
        tree as described for the corresponding pure virtual function in the base class. For this
        class, TO DO. */
    vector<TreeNode*> constructTree(TreeNode* root) override;
};

//////////////////////////////////////////////////////////////////////

#endif

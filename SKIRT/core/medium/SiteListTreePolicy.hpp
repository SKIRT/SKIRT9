/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SITELISTTREEPOLICY_HPP
#define SITELISTTREEPOLICY_HPP

#include "TreePolicy.hpp"
class Medium;
class Random;

//////////////////////////////////////////////////////////////////////

/** SiteListTreePolicy represents the configurable options and the corresponding implementation
    mechanisms for constructing spatial tree grids based on the positions of particles, sites or
    cells in an imported medium distribution. The policy locates a site list offered by one of the
    media components in the medium system. See the ImportedMedium class and the various Snapshot
    subclasses for more information on the site positions returned by imported media. If multiple
    media offer a site list, the first one in configuration order (i.e. ski file order) is used.

    In a first step the tree is subdivided in such a way that each leaf node contains at most one
    of the sites in the list. Subsequently each of these leaf nodes is further subdivided a fixed
    number of times, as configured by the user. However, the minimum and maximum tree subdvision
    levels (actually offered by the base class) override the other subdvision criteria described
    above. Tree nodes are always subdivided up to the minimum level, and nodes are never subdivided
    beyond the maximum level. */
class SiteListTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(SiteListTreePolicy, TreePolicy,
                  "a tree grid construction policy using positions defined by an imported medium")
        ATTRIBUTE_TYPE_ALLOWED_IF(SiteListTreePolicy, "SiteListInterface")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SiteListTreePolicy, "Level2")

        PROPERTY_INT(numExtraLevels, "the number of additional subdivision levels")
        ATTRIBUTE_MIN_VALUE(numExtraLevels, "0")
        ATTRIBUTE_MAX_VALUE(numExtraLevels, "30")
        ATTRIBUTE_DEFAULT_VALUE(numExtraLevels, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function constructs the hierarchical tree and all (interconnected) nodes forming the
        tree as described for the corresponding pure virtual function in the base class. The
        implementation for this class proceeds as decribed in the class header. After locating an
        object of a Medium subclass that offers the SiteListInterface, the tree is subdivided so
        that each leaf node contains at most one of the sites in the list. Subsequently each of
        these leaf nodes is further subdivided a fixed number of times, as indicated by the \em
        numExtraLevels attribute. */
    vector<TreeNode*> constructTree(TreeNode* root) override;
};

//////////////////////////////////////////////////////////////////////

#endif

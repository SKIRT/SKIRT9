/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREESPATIALGRIDTOPOLOGYPROBE_HPP
#define TREESPATIALGRIDTOPOLOGYPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** TreeSpatialGridTopologyProbe outputs a text data file representing the topology (i.e. the
    subdivision structure) of the tree grid used to discretize the spatial domain in the
    simulation. The resulting data file can be loaded by a FileTreeSpatialGrid instance to recreate
    a grid in another simulation that is identical to the one in the current simulation. This can
    be useful in situations where multiple simulations are being performed on input models with the
    same (or a very similar) spatial distribution of the transfer medium. Constructing a tree grid
    based on the medium density distribution takes time (because of the large number of density
    samples required), and the tree is likely to differ slightly between various runs (because the
    density is sampled at random positions, and the random number sequence varies because of
    parallelization). Loading the tree from the topology data file is much faster and guarantees
    that all simulations use identical tree grids. Note that the topology data stored in the file
    is scale-free; the simulation loading the topology must define the extent of the spatial
    domain.

    The output file is called <tt>prefix_probe_treetop.dat</tt>. After a brief descriptive header,
    the file contains lines with just a single integer number. The first line specifies the number
    of children for each nonleaf node (2 for a binary tree, 8 for an octtree, or 0 if the root node
    is not subdivided). The second line contains 1 if the root node is subdivided, or 0 if not. The
    following lines similarly contain 1 or 0 indicating subdivision for any children of the
    preceding node, recursively, in a depth-first traversal of the tree. */
class TreeSpatialGridTopologyProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(TreeSpatialGridTopologyProbe, SpecialtyProbe,
                  "specialty: data file representing the topology of the tree spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TreeSpatialGridTopologyProbe, "Level3&TreeSpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

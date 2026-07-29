/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHERICALCLUMPBVH_HPP
#define SPHERICALCLUMPBVH_HPP

#include "Box.hpp"
#include "Position.hpp"

////////////////////////////////////////////////////////////////////

/** SphericalClumpBVH is a low-level utility class that organizes a set of non-overlapping spheres
    ("clumps"), each defined by a center and a radius, in a data structure that allows efficient
    queries: finding the clump (if any) containing a given position, finding the nearest clump
    intersected by a ray, and so on.

    The implementation employs a linearized Bounding Volume Hierarchy (BVH) that is bulk-loaded
    using a surface area heuristic (SAH) for optimal balance. Construction of the data structure
    runs in a single serial thread. Because the data structure does not change after the initial
    call to loadClumps(), all query functions are thread-safe. */
class SphericalClumpBVH
{
public:
    /** An instance of this class represents a single spherical clump, defined by its center and
        radius. */
    class Clump
    {
        double _x, _y, _z, _r;

    public:
        /** This constructor initializes the clump to the given center coordinates and radius. */
        Clump(double x, double y, double z, double r) : _x{x}, _y{y}, _z{z}, _r{r} {}

        /** This function returns the center of the clump. */
        Position center() const { return Position(_x, _y, _z); }

        /** This function returns the radius of the clump. */
        double radius() const { return _r; }

        /** This function returns the Cartesian bounding box of the clump. */
        Box bounds() const { return Box(_x - _r, _y - _r, _z - _r, _x + _r, _y + _r, _z + _r); }
    };

    /** This function bulk-loads the specified clumps into the BVH, replacing any previous content.
        The clumps are not required to be mutually disjoint; see allDisjointClumps(). The specified
        vector must remain valid, unchanged, and at the same memory location for as long as the BVH
        is used, because the BVH stores a pointer to it rather than a copy. */
    void loadClumps(const vector<Clump>& clumps);

    /** This function returns the indices, into the vector most recently passed to loadClumps(), of
        a maximal subset of mutually non-overlapping clumps. The subset is built greedily in order
        of increasing index: clump 0 is always kept, and each subsequent clump is kept unless it
        overlaps a clump already kept. */
    vector<int> allDisjointClumps() const;

    /** This function returns the index, into the vector most recently passed to loadClumps(), of
        any clump containing the given position, or -1 if none does. If the loaded clumps are not
        mutually disjoint, which of the overlapping clumps is returned is unspecified. */
    int anyClumpContaining(Vec bfr) const;

    /** This function returns the index, into the vector most recently passed to loadClumps(), of
        the nearest clump intersected by the ray \f$({\bf{r}}_0,{\bf{k}})\f$ at a forward distance
        strictly between 0 and the value of sBest on entry, or -1 if there is none. On success, the
        distance to the entry point of that clump is stored in sBest. */
    int nearestClumpAlongRay(Position r0, Direction k, double& sBest) const;

    /** This function returns, for each clump that actually crosses the coordinate plane where
        coordinate \em axis (0 for x, 1 for y, 2 for z) equals \em value, a pair holding the index
        of that clump (into the vector most recently passed to loadClumps()) and the radius of the
        circle formed by the intersection of the clump's sphere with the plane. A cheap bounding-box
        test is used to prune the search, but every returned clump has already been confirmed by an
        exact sphere/plane test, so the caller does not need to repeat it. The ordering of the
        returned clumps is unspecified. */
    vector<std::pair<int, double>> clumpsCrossingPlane(int axis, double value) const;

    /** An instance of this struct holds summary statistics gathered by leafStatistics() from a
        walk of the leaves of an already-built BVH: for each of the leaf depth, the number of
        clumps in the leaf, and the diagonal of the leaf's bounding box, the minimum, maximum and
        average value taken over all leaves. */
    struct LeafStatistics
    {
        int minDepth{0};
        int maxDepth{0};
        double avgDepth{0.};
        int minCount{0};
        int maxCount{0};
        double avgCount{0.};
        double minDiagonal{0.};
        double maxDiagonal{0.};
        double avgDiagonal{0.};
    };

    /** This function walks all leaves of the already-built BVH and returns summary statistics on
        (1) the depth of each leaf in the tree (the root has depth zero), (2) the number of clumps
        held by each leaf, and (3) the diagonal of each leaf's bounding box. This is a read-only,
        optional diagnostic pass performed after construction -- it does not affect, and is not
        affected by, the build itself -- intended for logging or for tuning the BVH build
        parameters (NumBins and MaxLeafSize in the source file). If the BVH is empty (no clumps
        loaded), all statistics are zero. */
    LeafStatistics leafStatistics() const;

private:
    class Node
    {
    public:
        Box box;         // union of the bounding boxes of all clumps in this node's subtree
        int left;        // index into _nodes of the left child, or -1 if this is a leaf
        int right;       // index into _nodes of the right child, or -1 if this is a leaf
        int firstIndex;  // for a leaf: start of this leaf's range in _index
        int numIndices;  // for a leaf: number of clumps in this leaf; 0 for an interior node

        Node(const Box& box_, int left_, int right_, int firstIndex_, int numIndices_)
            : box(box_), left(left_), right(right_), firstIndex(firstIndex_), numIndices(numIndices_)
        {}

        bool isLeaf() const { return numIndices > 0; }
    };

    // flat array of BVH nodes; index 0 is the root, meaningful only when _numEntities > 0
    vector<Node> _nodes;

    // entity indices, reordered during the build so that the clumps of any single leaf
    // occupy the contiguous range [firstIndex, firstIndex+numIndices) of this array
    vector<int> _index;

    // pointer to array with clumps, passed on when bulk-loading the BVH
    const vector<Clump>* _clumps;

    // Recursively builds the BVH using binned SAH over the clumps referenced by
    // _index[begin,end), reordering that sub-range in place. Appends the newly created
    // node(s) to _nodes and returns the index of the node representing this sub-range.
    int buildRecursive(int begin, int end);
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SphericalClumpBVH.hpp"
#include "Quadrics.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // Number of candidate split planes ("bins") evaluated per axis during the SAH search
    // (binned SAH, Wald & Havran); O(N) work per node rather than O(N log N) for an exact
    // search, at negligible cost to tree quality.
    constexpr int NumBins = 256;

    // A leaf never holds more clumps than this, regardless of what SAH suggests. Bounds the
    // cost of the linear scan at each leaf, and guarantees the build terminates even when
    // many clumps have coincident or near-coincident centers.
    constexpr int MaxLeafSize = 32;
}

////////////////////////////////////////////////////////////////////

int SphericalClumpBVH::buildRecursive(int begin, int end)
{
    // aggregate bounding box for this sub-range: the union of the entity boxes referenced
    // by _index[begin,end).
    Box box = (*_clumps)[_index[begin]].bounds();
    for (int i = begin + 1; i != end; ++i) box.extend((*_clumps)[_index[i]].bounds());

    int numPrims = end - begin;
    int nodeIndex = static_cast<int>(_nodes.size());

    auto makeLeaf = [&]() {
        _nodes.emplace_back(box, -1, -1, begin, numPrims);
        return nodeIndex;
    };

    if (numPrims <= MaxLeafSize) return makeLeaf();

    // centroid bounding box, used to pick a split axis and the bin boundaries along it
    Position c0 = (*_clumps)[_index[begin]].center();
    double cmin[3] = {c0.x(), c0.y(), c0.z()};
    double cmax[3] = {c0.x(), c0.y(), c0.z()};
    for (int i = begin + 1; i != end; ++i)
    {
        Position c = (*_clumps)[_index[i]].center();
        double v[3] = {c.x(), c.y(), c.z()};
        for (int a = 0; a != 3; ++a)
        {
            cmin[a] = min(cmin[a], v[a]);
            cmax[a] = max(cmax[a], v[a]);
        }
    }
    int axis = 0;
    double extentOnAxis = cmax[0] - cmin[0];
    for (int a = 1; a != 3; ++a)
    {
        double e = cmax[a] - cmin[a];
        if (e > extentOnAxis)
        {
            extentOnAxis = e;
            axis = a;
        }
    }

    // all centroids coincide -> no split can separate the clumps; make a leaf
    // (also protects the binning step below from a division by zero)
    if (extentOnAxis <= 0.0) return makeLeaf();

    auto centroidAxis = [axis](Position c) { return axis == 0 ? c.x() : axis == 1 ? c.y() : c.z(); };

    // bin the clumps of [begin,end) along the chosen axis
    struct Bin
    {
        int count = 0;
        Box box;
        bool hasBox = false;
    };
    vector<Bin> bins(NumBins);
    double lo = cmin[axis], hi = cmax[axis];
    auto binIndexOf = [&](double v) {
        int b = static_cast<int>(NumBins * (v - lo) / (hi - lo));
        return min(max(b, 0), NumBins - 1);
    };
    for (int i = begin; i != end; ++i)
    {
        const Clump& c = (*_clumps)[_index[i]];
        int b = binIndexOf(centroidAxis(c.center()));
        Box eb = c.bounds();
        if (bins[b].hasBox)
            bins[b].box.extend(eb);
        else
        {
            bins[b].box = eb;
            bins[b].hasBox = true;
        }
        bins[b].count++;
    }

    // running unions/counts from the left and from the right, used to evaluate the SAH
    // cost of splitting after each bin boundary without re-scanning for every candidate
    vector<Box> prefixBox(NumBins);
    vector<int> prefixCount(NumBins);
    vector<Box> suffixBox(NumBins);
    vector<int> suffixCount(NumBins);
    {
        bool has = false;
        Box running;
        int runningCount = 0;
        for (int b = 0; b != NumBins; ++b)
        {
            if (bins[b].hasBox)
            {
                if (has)
                    running.extend(bins[b].box);
                else
                {
                    running = bins[b].box;
                    has = true;
                }
            }
            runningCount += bins[b].count;
            prefixBox[b] = running;
            prefixCount[b] = runningCount;
        }
    }
    {
        bool has = false;
        Box running;
        int runningCount = 0;
        for (int b = NumBins - 1; b >= 0; --b)
        {
            if (bins[b].hasBox)
            {
                if (has)
                    running.extend(bins[b].box);
                else
                {
                    running = bins[b].box;
                    has = true;
                }
            }
            runningCount += bins[b].count;
            suffixBox[b] = running;
            suffixCount[b] = runningCount;
        }
    }

    // evaluate the SAH cost of splitting right after bin b, for every interior bin boundary,
    // and keep the cheapest; cost is expressed relative to the parent's surface area so it
    // can be compared directly against the cost of not splitting at all
    double parentArea = box.surfaceArea();
    double bestCost = std::numeric_limits<double>::infinity();
    int bestSplit = -1;
    for (int b = 0; b != NumBins - 1; ++b)
    {
        if (prefixCount[b] == 0 || suffixCount[b + 1] == 0) continue;
        double cost =
            (prefixCount[b] * prefixBox[b].surfaceArea() + suffixCount[b + 1] * suffixBox[b + 1].surfaceArea())
            / parentArea;
        if (cost < bestCost)
        {
            bestCost = cost;
            bestSplit = b;
        }
    }

    // no viable split, or splitting is not cheaper than a leaf -> make a leaf
    double leafCost = static_cast<double>(numPrims);
    if (bestSplit < 0 || bestCost >= leafCost) return makeLeaf();

    double splitPos = lo + (hi - lo) * (bestSplit + 1) / NumBins;
    auto midIt = std::partition(_index.begin() + begin, _index.begin() + end,
                                [&](int m) { return centroidAxis((*_clumps)[m].center()) < splitPos; });
    int mid = static_cast<int>(midIt - _index.begin());
    if (mid == begin || mid == end) mid = (begin + end) / 2;  // guard against a degenerate partition

    // reserve this node's slot now; children are appended afterwards and get higher
    // indices, so we come back and fill in left/right once the recursion returns
    _nodes.emplace_back(box, -1, -1, 0, 0);
    int leftIndex = buildRecursive(begin, mid);
    int rightIndex = buildRecursive(mid, end);
    _nodes[nodeIndex].left = leftIndex;
    _nodes[nodeIndex].right = rightIndex;
    return nodeIndex;
}

////////////////////////////////////////////////////////////////////

void SphericalClumpBVH::loadClumps(const vector<Clump>& clumps)
{
    _clumps = &clumps;
    int numClumps = static_cast<int>(clumps.size());

    _nodes.clear();
    _index.resize(numClumps);
    for (int m = 0; m != numClumps; ++m) _index[m] = m;

    if (numClumps > 0)
    {
        // rough heuristic upper bound on the final node count (a binary tree with roughly
        // numEntities/MaxLeafSize leaves has on the order of twice as many nodes total)
        _nodes.reserve(2 * numClumps / MaxLeafSize + 1);
        buildRecursive(0, numClumps);
    }
}

////////////////////////////////////////////////////////////////////

vector<int> SphericalClumpBVH::allDisjointClumps() const
{
    vector<int> result;
    int numClumps = static_cast<int>(_index.size());
    if (numClumps == 0) return result;

    // accepted[m] becomes true once entity m has been added to the result; a plain
    // vector<char> is used instead of vector<bool> for fast, unambiguous random-access.
    vector<char> accepted(numClumps, 0);
    accepted[0] = 1;
    result.push_back(0);

    // For each remaining clump, we need to know whether it overlaps ANY already-accepted clump
    vector<int> stack;
    for (int i = 1; i != numClumps; ++i)
    {
        const Clump& ci = (*_clumps)[i];
        Box queryBox = ci.bounds();
        bool conflict = false;

        stack.clear();
        stack.push_back(0);
        while (!stack.empty() && !conflict)
        {
            int nodeIndex = stack.back();
            stack.pop_back();
            const Node& node = _nodes[nodeIndex];
            if (!node.box.intersects(queryBox)) continue;

            if (node.isLeaf())
            {
                for (int k = node.firstIndex; k != node.firstIndex + node.numIndices; ++k)
                {
                    int m = _index[k];
                    if (!accepted[m]) continue;
                    const Clump& cm = (*_clumps)[m];
                    if (Quadrics::doSpheresOverlap(ci.center(), ci.radius(), cm.center(), cm.radius()))
                    {
                        conflict = true;
                        break;
                    }
                }
            }
            else
            {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }

        if (!conflict)
        {
            accepted[i] = 1;
            result.push_back(i);
        }
    }
    return result;
}

////////////////////////////////////////////////////////////////////

int SphericalClumpBVH::anyClumpContaining(Vec bfr) const
{
    if (_nodes.empty()) return -1;

    // use a thread_local stack to avoid allocations between consecutive queries
    thread_local vector<int> stack;
    stack.clear();
    stack.push_back(0);

    while (!stack.empty())
    {
        int nodeIndex = stack.back();
        stack.pop_back();
        const Node& node = _nodes[nodeIndex];
        if (!node.box.contains(bfr)) continue;

        if (node.isLeaf())
        {
            for (int k = node.firstIndex; k != node.firstIndex + node.numIndices; ++k)
            {
                int m = _index[k];
                const Clump& c = (*_clumps)[m];
                if (Quadrics::isPositionInSphere(bfr, c.center(), c.radius())) return m;
            }
        }
        else
        {
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////

// This is a nearest-first, depth-first traversal with branch-and-bound pruning: a node's
// ray-box entry distance (from Box::intersects) is a valid lower bound on the entry
// distance of any clump inside it, so a node whose entry distance is not smaller than the
// best confirmed hit so far can be skipped outright. Visiting the geometrically nearer
// child first tends to tighten that bound quickly. Unlike an ordered multi-hit iterator,
// this only ever needs the single nearest hit, since the caller typically re-queries fresh
// from its new position after every hit.
int SphericalClumpBVH::nearestClumpAlongRay(Position r0, Direction k, double& sBest) const
{
    if (_nodes.empty()) return -1;

    int bestIndex = -1;

    struct StackEntry
    {
        int nodeIndex;
        double boxEntry;
    };
    thread_local vector<StackEntry> stack;  // thread_local to avoid reallocations on repeated use
    stack.clear();

    double smin, smax;
    if (_nodes[0].box.intersects(r0, k, smin, smax) && smin < sBest) stack.push_back({0, smin});

    while (!stack.empty())
    {
        StackEntry top = stack.back();
        stack.pop_back();
        if (top.boxEntry >= sBest) continue;  // superseded by a tighter bound found meanwhile

        const Node& node = _nodes[top.nodeIndex];
        if (node.isLeaf())
        {
            for (int i = node.firstIndex; i != node.firstIndex + node.numIndices; ++i)
            {
                int m = _index[i];
                const Clump& c = (*_clumps)[m];
                double s0 = Quadrics::firstIntersectionSphere(r0, k, c.center(), c.radius());
                if (s0 > 0. && s0 < sBest)
                {
                    sBest = s0;
                    bestIndex = m;
                }
            }
        }
        else
        {
            double sminL, smaxL, sminR, smaxR;
            bool hitL = _nodes[node.left].box.intersects(r0, k, sminL, smaxL) && sminL < sBest;
            bool hitR = _nodes[node.right].box.intersects(r0, k, sminR, smaxR) && sminR < sBest;

            // push the farther child first so the nearer one is popped (and processed) first
            if (hitL && hitR)
            {
                if (sminL <= sminR)
                {
                    stack.push_back({node.right, sminR});
                    stack.push_back({node.left, sminL});
                }
                else
                {
                    stack.push_back({node.left, sminL});
                    stack.push_back({node.right, sminR});
                }
            }
            else if (hitL)
                stack.push_back({node.left, sminL});
            else if (hitR)
                stack.push_back({node.right, sminR});
        }
    }

    return bestIndex;
}

////////////////////////////////////////////////////////////////////

vector<std::pair<int, double>> SphericalClumpBVH::clumpsCrossingPlane(int axis, double value) const
{
    vector<std::pair<int, double>> result;
    if (_nodes.empty()) return result;

    auto straddles = [axis, value](const Box& b) {
        double lo = axis == 0 ? b.xmin() : axis == 1 ? b.ymin() : b.zmin();
        double hi = axis == 0 ? b.xmax() : axis == 1 ? b.ymax() : b.zmax();
        return lo <= value && value <= hi;
    };

    vector<int> stack;
    stack.push_back(0);
    while (!stack.empty())
    {
        int nodeIndex = stack.back();
        stack.pop_back();
        const Node& node = _nodes[nodeIndex];
        if (!straddles(node.box)) continue;

        if (node.isLeaf())
        {
            for (int k = node.firstIndex; k != node.firstIndex + node.numIndices; ++k)
            {
                int m = _index[k];
                const Clump& c = (*_clumps)[m];
                Position center = c.center();
                double outOfPlaneCoord = (axis == 0 ? center.x() : axis == 1 ? center.y() : center.z()) - value;
                double circleRadius;
                if (Quadrics::sphereIntersectsPlane(c.radius(), outOfPlaneCoord, circleRadius))
                    result.emplace_back(m, circleRadius);
            }
        }
        else
        {
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }
    return result;
}

////////////////////////////////////////////////////////////////////

SphericalClumpBVH::LeafStatistics SphericalClumpBVH::leafStatistics() const
{
    LeafStatistics stats;
    if (_nodes.empty()) return stats;

    stats.minDepth = std::numeric_limits<int>::max();
    stats.minCount = std::numeric_limits<int>::max();
    stats.minDiagonal = std::numeric_limits<double>::infinity();

    int numLeaves = 0;
    long depthSum = 0;
    long countSum = 0;
    double diagonalSum = 0.;

    struct StackEntry
    {
        int nodeIndex;
        int depth;
    };
    vector<StackEntry> stack;
    stack.push_back({0, 0});
    while (!stack.empty())
    {
        StackEntry top = stack.back();
        stack.pop_back();
        const Node& node = _nodes[top.nodeIndex];

        if (node.isLeaf())
        {
            numLeaves++;
            depthSum += top.depth;
            countSum += node.numIndices;
            double diagonal = node.box.diagonal();
            diagonalSum += diagonal;

            stats.minDepth = min(stats.minDepth, top.depth);
            stats.maxDepth = max(stats.maxDepth, top.depth);
            stats.minCount = min(stats.minCount, node.numIndices);
            stats.maxCount = max(stats.maxCount, node.numIndices);
            stats.minDiagonal = min(stats.minDiagonal, diagonal);
            stats.maxDiagonal = max(stats.maxDiagonal, diagonal);
        }
        else
        {
            stack.push_back({node.left, top.depth + 1});
            stack.push_back({node.right, top.depth + 1});
        }
    }

    stats.avgDepth = static_cast<double>(depthSum) / numLeaves;
    stats.avgCount = static_cast<double>(countSum) / numLeaves;
    stats.avgDiagonal = diagonalSum / numLeaves;
    return stats;
}

////////////////////////////////////////////////////////////////////

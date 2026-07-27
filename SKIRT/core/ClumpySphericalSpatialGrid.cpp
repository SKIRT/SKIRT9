/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClumpySphericalSpatialGrid.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Quadrics.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

// Geometric helper functions
namespace
{
    // small value, used to progress a path; leave at 1e-10 -- with 1e-11, inaccuracies start to occur
    constexpr double EPS = 1e-10;

    // writes the axisymmetric structured-grid lines (spheres and cones) to a meridional plot;
    // shared by write_xz() and write_yz() since the structured grid itself is axisymmetric
    // (only the overlaid clumps differ between those two views)
    void writeMeridionalStructure(SpatialGridPlotFile* outfile, const Array& rv, const Array& thetav)
    {
        for (double r : rv) outfile->writeCircle(r);

        double r0 = rv[0];
        double r1 = rv[rv.size() - 1];
        for (double theta : thetav)
        {
            outfile->writeLine(r0 * sin(theta), r0 * cos(theta), r1 * sin(theta), r1 * cos(theta));
            outfile->writeLine(-r0 * sin(theta), -r0 * cos(theta), -r1 * sin(theta), -r1 * cos(theta));
        }
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // Number of candidate split planes ("bins") evaluated per axis during the SAH search
    // (binned SAH, Wald & Havran); O(N) work per node rather than O(N log N) for an exact
    // search, at negligible cost to tree quality.
    constexpr int NumBins = 16;

    // A leaf never holds more clumps than this, regardless of what SAH suggests. Bounds the
    // cost of the linear scan at each leaf, and guarantees the build terminates even when
    // many clumps have coincident or near-coincident centers.
    constexpr int MaxLeafSize = 4;
}

////////////////////////////////////////////////////////////////////

// Private utility class for organizing spheres in a data structure that allows efficient queries.
// The implementation employs a linearized Bounding Volume Hierarchy (BVH) that is bulk-loaded
// using a surface area heuristic (SAH) for optimal balance. Construction of the data structure
// runs in a single serial thread. Because the data structure does not change after initial
// construction, all queries are thread-safe.
class ClumpySphericalSpatialGrid::BVH
{
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
    int buildRecursive(int begin, int end)
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

public:
    // Bulk-loads the specified clumps into the BVH
    void loadClumps(const vector<Clump>& clumps)
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

    // Returns the indices of a maximal subset of mutually non-overlapping clumps, built
    // greedily in order of increasing index: clump 0 is always kept, and each subsequent
    // clump is kept unless it overlaps a clump already kept. Uses the BVH already built
    // over ALL clumps (including ones that will be rejected here) purely to prune the
    // overlap search down to nearby candidates via box intersection, rather than checking
    // against every previously-accepted clump directly (which would be O(M^2) overall).
    vector<int> allDisjointClumps() const
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

    // Returns the index of any clump containing the given position, or -1 if none does.
    // No ordering is needed (clumps are non-overlapping by the time this is used in
    // practice), so a plain box-pruned depth-first search returning on the first hit is
    // the fastest option.
    int anyClumpContaining(Vec bfr) const
    {
        if (_nodes.empty()) return -1;

        // use a thread_local stack to avoid allocations between consecutibe queries
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

    // Returns the index of the nearest clump intersected by the ray (r0,k) at a forward
    // distance strictly between 0 and sBest, or -1 if there is none. On success, sBest is
    // set to the entry distance of that specific clump.
    //
    // This is a nearest-first, depth-first traversal with branch-and-bound pruning: a node's
    // ray-box entry distance (from Box::intersects) is a valid lower bound on the entry
    // distance of any clump inside it, so a node whose entry distance is not smaller than the
    // best confirmed hit so far can be skipped outright. Visiting the geometrically nearer
    // child first tends to tighten that bound quickly. Unlike an ordered multi-hit iterator,
    // this only ever needs the single nearest hit, since the segment generator re-queries
    // fresh from its new position after every segment.
    int nearestClumpAlongRay(Position r0, Direction k, double& sBest) const
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

    // Returns the indices of the clumps whose bounding box may cross the coordinate plane
    // where coordinate 'axis' (0=x, 1=y, 2=z) equals 'value'. Runs once at plotting time, so
    // a plain box-pruned traversal (no thread_local reuse) is fine. The result can include a
    // few clumps whose bounding box straddles the plane but whose actual sphere doesn't quite
    // reach it; the caller is expected to apply the exact test itself, which it needs anyway
    // to get the intersection circle's radius.
    vector<int> clumpsCrossingPlane(int axis, double value) const
    {
        vector<int> result;
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
                for (int k = node.firstIndex; k != node.firstIndex + node.numIndices; ++k) result.push_back(_index[k]);
            }
            else
            {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }
        return result;
    }
};

////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::setupSelfAfter()
{
    SphereSpatialGrid::setupSelfAfter();
    Log* log = find<Log>();

    // ---- build the structured grid (copied from Sphere3DSpatialGrid::setupSelfAfter) ----

    // radial
    _Nr = _meshRadial->numBins();
    _rv = _meshRadial->mesh() * (maxRadius() - minRadius()) + minRadius();

    // limit the epsilon we use for progressing the path to a value smaller than the hole and/or the first bin
    bool hasHole = _rv[0] > 0.;
    double meshEps = min(EPS * maxRadius(), 0.1 * (hasHole ? min(_rv[0], _rv[1] - _rv[0]) : _rv[1]));

    // if the grid has no hole, create an artificial hole larger than the epsilon we use for progressing the path
    // so that the segment generator has a chance to reset the phi bin index when crossing the origin
    if (!hasHole) _rv[0] = 2. * meshEps;

    // polar
    _Ntheta = initPolarGrid(_thetav, _cv, _meshPolar);

    // azimuthal
    _Nphi = _meshAzimuthal->numBins();
    _phiv = _meshAzimuthal->mesh() * (2. * M_PI) - M_PI;

    // verify that the azimuth bins are smaller than 120 degrees, with some leeway,
    // so that the segment generator never inadvertently intersects the path with the reflected phi bin border
    constexpr double maxPhi = 0.7 * M_PI;
    for (int k = 0; k != _Nphi; ++k)
        if (_phiv[k + 1] - _phiv[k] > maxPhi) throw FATALERROR("Azimuth bin is wider than 120 deg");

    // pre-calculate sines and cosines for azimuthal bin borders; make sure that the outer boundary values are exact
    _sinv = sin(_phiv);
    _cosv = cos(_phiv);
    _sinv[0] = 0.;
    _cosv[0] = -1.;
    _sinv[_Nphi] = 0.;
    _cosv[_Nphi] = -1.;

    // number of structured cells
    _Ncells = _Nr * _Ntheta * _Nphi;

    // precompute the nominal (clump-free) volume of each structured cell
    _cellVolume.resize(_Ncells);
    for (int s = 0; s != _Ncells; ++s)
    {
        double rmin, thetamin, phimin, rmax, thetamax, phimax;
        getCoords(s, rmin, thetamin, phimin, rmax, thetamax, phimax);
        _cellVolume[s] = (1. / 3.) * Quadrics::pow3(rmin, rmax) * (cos(thetamin) - cos(thetamax)) * (phimax - phimin);
    }

    // ---- read the file defining the spherical clumps, rejecting those outside the domain ----

    int numOutside = 0;
    {
        TextInFile infile(this, _filename, "clump centers and radii");
        infile.addColumn("position x", "length", "pc");
        infile.addColumn("position y", "length", "pc");
        infile.addColumn("position z", "length", "pc");
        infile.addColumn("radius r", "length", "pc");
        double x, y, z, r;
        while (infile.readRow(x, y, z, r))
        {
            if (Quadrics::isSphereInShell(Vec(x, y, z), r, minRadius(), maxRadius()))
                _clumps.emplace_back(x, y, z, r);
            else
                numOutside++;
        }
        infile.close();
    }

    // build a BVH over all clumps inside the domain, including mutually overlapping ones
    _bvh = new BVH;
    log->info("Constructing bounding volume hierarchy for " + std::to_string(_clumps.size()) + " clumps...");
    _bvh->loadClumps(_clumps);

    // check for overlapping clumps using the BVH
    log->info("Checking for overlapping clumps...");
    vector<int> disjointIndices = _bvh->allDisjointClumps();
    int numOverlapping = _clumps.size() - disjointIndices.size();

    // if there are any overlapping clumps, rebuild the list of clumps and rebuild the BVH
    if (numOverlapping > 0)
    {
        log->info("Rebuilding bounding volume hierarchy for non-overlapping clumps...");
        vector<Clump> originalClumps;
        originalClumps.swap(_clumps);
        for (int m : disjointIndices) _clumps.emplace_back(originalClumps[m]);
        _bvh->loadClumps(_clumps);
    }

    // remember the final number of clumps, which also serves as the index offset of the structured cells
    _numClumps = _clumps.size();

    // ---- for each clump, Monte Carlo sample its overlap with the structured cells ----
    //
    // We sample points uniformly inside the clump's own sphere (rather than inside each candidate
    // structured cell) and determine which structured cell each sample lands in. Because every
    // sample lands in exactly one cell, the resulting per-cell hit counts always sum exactly to
    // the total number of samples drawn, so the fractional volume contributions derived from them
    // always sum exactly to the clump's true volume, with no leakage even though the individual
    // per-cell splits are Monte Carlo estimates.
    //
    // A single fixed sample budget per clump would grow increasingly inaccurate for a clump that
    // straddles many small cells, since each individual cell would then receive only a small
    // fraction of the samples. To avoid that, sampling proceeds in two phases: phase 1 draws K
    // samples and discovers how many distinct cells N are touched; if N>1, phase 2 tops up the
    // budget to a total of N*K samples, so that every touched cell receives on the order of K
    // samples on average, regardless of how many cells the clump straddles.

    Random* rnd = random();
    int K = find<Configuration>()->numDensitySamples();

    // draws numSamples additional samples uniformly inside the clump's sphere and tallies, for
    // each, which structured cell it lands in; typically only a handful of cells are touched, so
    // a linear-scan association list is lighter than a dense per-clump array of size _Ncells
    auto sampleAndTally = [&](const Clump& clump, int numSamples, vector<std::pair<int, int>>& hits) {
        for (int n = 0; n != numSamples; ++n)
        {
            Position p = rnd->positionInSphere(clump.center(), clump.radius());
            double r, theta, phi;
            p.spherical(r, theta, phi);
            int i = NR::locateClip(_rv, r);
            int j = NR::locateClip(_thetav, theta);
            int k = NR::locateClip(_phiv, phi);
            int s = structuredIndex(i, j, k);

            bool found = false;
            for (auto& hit : hits)
                if (hit.first == s)
                {
                    hit.second++;
                    found = true;
                    break;
                }
            if (!found) hits.emplace_back(s, 1);
        }
    };

    int numStraddling = 0;
    log->info("Estimating structured-cell volumes overlapped by " + std::to_string(_numClumps) + " clumps (at least "
              + std::to_string(K) + " samples per clump)...");
    vector<std::pair<int, int>> hits;  // reused across clumps to avoid reallocating on every iteration
    for (int ci = 0; ci != _numClumps; ++ci)
    {
        Clump& clump = _clumps[ci];
        double clumpVolume = Quadrics::volumeSphere(clump.radius());

        // phase 1: draw K samples and discover how many distinct cells are touched
        hits.clear();
        sampleAndTally(clump, K, hits);
        int totalSamples = K;

        // phase 2: if more than one cell was touched, top up the sample budget so that every
        // touched cell receives on the order of K samples on average
        if (hits.size() > 1)
        {
            int extraSamples = (static_cast<int>(hits.size()) - 1) * K;
            sampleAndTally(clump, extraSamples, hits);
            totalSamples += extraSamples;
        }

        int numOverlappingCells = static_cast<int>(hits.size());
        clump.setNumOverlappingCells(numOverlappingCells);
        if (numOverlappingCells > 1) numStraddling++;

        for (const auto& hit : hits)
        {
            double fraction = static_cast<double>(hit.second) / totalSamples;
            _cellVolume[hit.first] -= clumpVolume * fraction;
        }
    }

    // clamp away any small negative volumes caused by Monte Carlo noise or near-coincident boundaries
    for (double& v : _cellVolume) v = max(v, 0.);

    // inform the user
    log->info("Summary:");
    log->info("  Clumps read from file:      " + std::to_string(_numClumps + numOutside + numOverlapping));
    log->info("  Rejected (outside domain):  " + std::to_string(numOutside));
    log->info("  Rejected (overlapping):     " + std::to_string(numOverlapping));
    log->info("  Remaining in spatial grid:  " + std::to_string(_numClumps));
    log->info("  Straddling a structured cell wall: " + std::to_string(numStraddling));

    // ---- epsilon used to progress a path, small relative to both the mesh and the smallest clump ----

    double smallestClumpRadius = maxRadius();
    for (const Clump& clump : _clumps) smallestClumpRadius = min(smallestClumpRadius, clump.radius());
    double clumpEps = _numClumps ? min(EPS * maxRadius(), 0.1 * smallestClumpRadius) : meshEps;
    _eps = min(meshEps, clumpEps);
}

//////////////////////////////////////////////////////////////////////

ClumpySphericalSpatialGrid::~ClumpySphericalSpatialGrid()
{
    delete _bvh;
}

//////////////////////////////////////////////////////////////////////

int ClumpySphericalSpatialGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int ClumpySphericalSpatialGrid::numCells() const
{
    return _numClumps + _Ncells;
}

//////////////////////////////////////////////////////////////////////

double ClumpySphericalSpatialGrid::volume(int m) const
{
    if (m >= 0 && m < _numClumps) return Quadrics::volumeSphere(_clumps[m].radius());
    int s = m - _numClumps;
    if (s >= 0 && s < _Ncells) return _cellVolume[s];
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double ClumpySphericalSpatialGrid::diagonal(int m) const
{
    if (m >= 0 && m < _numClumps) return 2. * _clumps[m].radius();

    int s = m - _numClumps;
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(s, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        Position p0(rmin, thetamin, phimin, Position::CoordinateSystem::SPHERICAL);
        Position p1(rmax, thetamax, phimax, Position::CoordinateSystem::SPHERICAL);
        return (p1 - p0).norm();
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

int ClumpySphericalSpatialGrid::cellIndex(Position bfr) const
{
    int c = _bvh->anyClumpContaining(bfr);
    if (c >= 0) return c;

    double r, theta, phi;
    bfr.spherical(r, theta, phi);
    int i = NR::locateFail(_rv, r);
    if (i < 0) return -1;
    int j = NR::locateClip(_thetav, theta);
    int k = NR::locateClip(_phiv, phi);

    return _numClumps + structuredIndex(i, j, k);
}

//////////////////////////////////////////////////////////////////////

Position ClumpySphericalSpatialGrid::centralPositionInCell(int m) const
{
    if (m >= 0 && m < _numClumps) return _clumps[m].center();

    int s = m - _numClumps;
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(s, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        double r = 0.5 * (rmin + rmax);
        double theta = 0.5 * (thetamin + thetamax);
        double phi = 0.5 * (phimin + phimax);
        Position candidate(r, theta, phi, Position::CoordinateSystem::SPHERICAL);

        // if the nominal center happens to fall inside one of the clumps, fall back to a random
        // position instead; the BVH gives an exact answer, unlike the (possibly incomplete, since
        // it is built from Monte Carlo discovery -- see setupSelfAfter) per-cell clump list
        if (_bvh->anyClumpContaining(candidate) >= 0) return randomPositionInCell(m);
        return candidate;
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

Position ClumpySphericalSpatialGrid::randomPositionInCell(int m) const
{
    if (m >= 0 && m < _numClumps) return random()->positionInSphere(_clumps[m].center(), _clumps[m].radius());

    int s = m - _numClumps;
    double rmin, thetamin, phimin, rmax, thetamax, phimax;
    if (getCoords(s, rmin, thetamin, phimin, rmax, thetamax, phimax))
    {
        while (true)
        {
            double r = cbrt(Quadrics::pow3(rmin) + Quadrics::pow3(rmin, rmax) * random()->uniform());
            double theta = acos(cos(thetamin) + (cos(thetamax) - cos(thetamin)) * random()->uniform());
            double phi = phimin + (phimax - phimin) * random()->uniform();
            Position candidate(r, theta, phi, Position::CoordinateSystem::SPHERICAL);

            // accept the position unless it falls inside one of the clumps; the BVH gives an
            // exact answer, unlike the (possibly incomplete, since it is built from Monte Carlo
            // discovery -- see setupSelfAfter) per-cell clump list
            if (_bvh->anyClumpContaining(candidate) < 0) return candidate;
        }
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

class ClumpySphericalSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const ClumpySphericalSpatialGrid* _grid{nullptr};
    double _eps{0.};
    int _i{-1}, _j{-1}, _k{-1};  // structured bin indices; meaningful only while not inside a clump
    int _clump{-1};              // index of the clump currently inside, or -1 for the structured grid

public:
    MySegmentGenerator(const ClumpySphericalSpatialGrid* grid) : _grid(grid), _eps(grid->_eps) {}

    // determines and sets the indices i, j and k of the structured cell containing the current position
    //   i is set to -1 if the position is inside rmin and to Nr if the position is outside rmax
    //   j is clipped to the range 0..Ntheta-1
    //   k is clipped to the range 0..Nphi-1
    // returns true if the position is inside rmax, false if it is outside rmax
    bool setCellIndices()
    {
        double radius, theta, phi;
        r().spherical(radius, theta, phi);
        _i = NR::locate(_grid->_rv, radius);
        _j = NR::locateClip(_grid->_thetav, theta);
        _k = NR::locateClip(_grid->_phiv, phi);
        return _i < _grid->_Nr;
    }

    // sets the state to outside and returns false
    bool abortPath()
    {
        setState(State::Outside);
        return false;
    }

    // returns the distance to the first intersection between the current path and the radial
    // boundary sphere with given bin index, or 0 if there is none
    double firstIntersectionRadialBoundary(int i) { return Quadrics::firstIntersectionSphere(r(), k(), _grid->_rv[i]); }

    // returns the distance to the first intersection between the current path and the cone with
    // given bin index, or 0 if there is none (the degenerate cone with zero cosine is treated separately)
    double firstIntersectionCone(int j) { return Quadrics::firstIntersectionCone(r(), k(), _grid->_cv[j]); }

    // returns the distance to the intersection between the current path and the meridional plane
    // with the given bin index, or 0 if there is none
    double intersectionMeridionalPlane(int kbin)
    {
        return Quadrics::intersectionMeridionalPlane(r(), k(), _grid->_sinv[kbin], _grid->_cosv[kbin]);
    }

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // if necessary, try moving the path inside the grid
                if (r().norm() > _grid->maxRadius())
                {
                    double ds = firstIntersectionRadialBoundary(_grid->_Nr);
                    if (ds <= 0.) return abortPath();

                    propagater(ds + _eps);
                    if (!setCellIndices()) return abortPath();

                    setEmptySegment(ds);
                    setState(State::Inside);
                    return true;
                }

                // the original position was inside the grid; it may already be inside a clump
                // (e.g. a scattered photon relaunched from a scattering event inside a clump)
                _clump = _grid->_bvh->anyClumpContaining(r());
                if (!setCellIndices()) return abortPath();
                setState(State::Inside);

                // intentionally falls through to determine the first actual segment
            }

            case State::Inside:
            {
                // if currently inside a clump, the segment ends at the sphere boundary
                if (_clump >= 0)
                {
                    const Clump& c = _grid->_clumps[_clump];
                    double ds = Quadrics::firstIntersectionSphere(r(), k(), c.center(), c.radius());
                    if (ds <= 0.) return abortPath();
                    setSegment(_clump, ds);
                    propagater(ds + _eps);
                    _clump = -1;

                    // the clump may have carried the path across a structured-cell wall while
                    // inside it, so the cell indices cannot be trusted and must be recomputed
                    if (!setCellIndices()) return abortPath();
                    return true;
                }

                // if we're not inside the real or artificial hole, proceed to the next boundary in the regular way
                if (_i >= 0)
                {
                    // remember the indices of the current cell
                    int icur = _i;
                    int jcur = _j;
                    int kcur = _k;

                    // calculate the distance travelled inside the cell by considering the potential
                    // exit points for each of the six cell boundaries; the smallest positive
                    // intersection "distance" wins
                    double ds = DBL_MAX;

                    // inner radial boundary (always nonzero)
                    {
                        double s = firstIntersectionRadialBoundary(icur);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur - 1;  // may be decremented to -1 (inside the innermost boundary)
                            _j = jcur;
                            _k = kcur;
                        }
                    }

                    // outer radial boundary
                    {
                        double s = firstIntersectionRadialBoundary(icur + 1);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur + 1;  // may be incremented to Nr (beyond the outermost boundary)
                            _j = jcur;
                            _k = kcur;
                        }
                    }

                    // upper angular boundary (not applicable to uppermost cell)
                    if (jcur > 0)
                    {
                        double s = firstIntersectionCone(jcur);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur - 1;
                            _k = kcur;
                        }
                    }

                    // lower angular boundary (not applicable to lowest cell)
                    if (jcur < _grid->_Ntheta - 1)
                    {
                        double s = firstIntersectionCone(jcur + 1);
                        if (s > 0 && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur + 1;
                            _k = kcur;
                        }
                    }

                    // clockwise azimuthal boundary
                    {
                        double s = intersectionMeridionalPlane(kcur);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur;
                            _k = kcur > 0 ? kcur - 1 : _grid->_Nphi - 1;  // scroll from -pi to pi
                        }
                    }

                    // anticlockwise azimuthal boundary
                    {
                        double s = intersectionMeridionalPlane(kcur + 1);
                        if (s > 0. && s < ds)
                        {
                            ds = s;
                            _i = icur;
                            _j = jcur;
                            _k = (kcur + 1) % _grid->_Nphi;  // scroll from pi to -pi
                        }
                    }

                    // if no exit point was found, abort the path
                    if (ds == DBL_MAX) return abortPath();

                    // check whether a clump is entered before this structured-cell boundary is reached;
                    // bounding the BVH search by ds both prunes the search and means "no clump found"
                    // and "found beyond the cell boundary" collapse into the same outcome
                    double dsClump = ds;
                    int clumpHit = _grid->_bvh->nearestClumpAlongRay(r(), k(), dsClump);
                    if (clumpHit >= 0)
                    {
                        // the segment up to the clump entry belongs to the current structured cell;
                        // the tentative (i,j,k) computed above never actually apply since that
                        // boundary is not reached -- they are simply recomputed from scratch once
                        // the clump is exited, so there is no need to restore them here
                        setSegment(_grid->_numClumps + _grid->structuredIndex(icur, jcur, kcur), dsClump);
                        propagater(dsClump + _eps);
                        _clump = clumpHit;
                    }
                    else
                    {
                        setSegment(_grid->_numClumps + _grid->structuredIndex(icur, jcur, kcur), ds);
                        propagater(ds + _eps);
                        if (_i >= _grid->_Nr) setState(State::Outside);
                    }
                }

                // if we're inside the hole, skip to the hole radius in one empty segment step
                // and recalculate the bin indices (the phi bin index changes when crossing the origin)
                else
                {
                    double ds = firstIntersectionRadialBoundary(0);
                    if (ds <= 0.) return abortPath();
                    setEmptySegment(ds);
                    propagater(ds + _eps);
                    if (!setCellIndices()) return abortPath();
                }
                return true;
            }

            case State::Outside:
            {
            }
        }
        return false;  // unreachable; silences a compiler warning about a missing return
    }
};

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> ClumpySphericalSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

//////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // structured grid: spheres and meridional planes
    for (double r : _rv) outfile->writeCircle(r);
    for (double phi : _phiv)
        outfile->writeLine(_rv[0] * cos(phi), _rv[0] * sin(phi), _rv[_Nr] * cos(phi), _rv[_Nr] * sin(phi));

    // spherical clumps crossing the xy plane
    for (int m : _bvh->clumpsCrossingPlane(2, 0.))
    {
        const Clump& c = _clumps[m];
        double circleRadius;
        if (Quadrics::sphereIntersectsPlane(c.radius(), c.center().z(), circleRadius))
            outfile->writeCircle(c.center().x(), c.center().y(), circleRadius);
    }
}

//////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    // structured grid: spheres and cones
    writeMeridionalStructure(outfile, _rv, _thetav);

    // spherical clumps crossing the xz plane
    for (int m : _bvh->clumpsCrossingPlane(1, 0.))
    {
        const Clump& c = _clumps[m];
        double circleRadius;
        if (Quadrics::sphereIntersectsPlane(c.radius(), c.center().y(), circleRadius))
            outfile->writeCircle(c.center().x(), c.center().z(), circleRadius);
    }
}

//////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    // structured grid: spheres and cones (identical to the xz view since the structured grid
    // itself is axisymmetric; unlike Sphere3DSpatialGrid we cannot simply delegate to write_xz()
    // for the full view here, because the overlaid clumps are not azimuthally symmetric)
    writeMeridionalStructure(outfile, _rv, _thetav);

    // spherical clumps crossing the yz plane
    for (int m : _bvh->clumpsCrossingPlane(0, 0.))
    {
        const Clump& c = _clumps[m];
        double circleRadius;
        if (Quadrics::sphereIntersectsPlane(c.radius(), c.center().x(), circleRadius))
            outfile->writeCircle(c.center().y(), c.center().z(), circleRadius);
    }
}

//////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // intentionally unimplemented
    (void)outfile;
}

//////////////////////////////////////////////////////////////////////

int ClumpySphericalSpatialGrid::structuredIndex(int i, int j, int k) const
{
    return k + (j + i * _Ntheta) * _Nphi;
}

//////////////////////////////////////////////////////////////////////

bool ClumpySphericalSpatialGrid::getCoords(int s, double& rmin, double& thetamin, double& phimin, double& rmax,
                                           double& thetamax, double& phimax) const
{
    if (s < 0 || s >= _Ncells) return false;

    int i = s / (_Ntheta * _Nphi);
    int j = (s / _Nphi) % _Ntheta;
    int k = s % _Nphi;

    rmin = _rv[i];
    thetamin = _thetav[j];
    phimin = _phiv[k];
    rmax = _rv[i + 1];
    thetamax = _thetav[j + 1];
    phimax = _phiv[k + 1];
    return true;
}

//////////////////////////////////////////////////////////////////////

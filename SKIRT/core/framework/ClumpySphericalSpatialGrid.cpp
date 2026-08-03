/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClumpySphericalSpatialGrid.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "PathSegmentGenerator.hpp"
#include "Quadrics.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of clumps drawn in the 3D grid plot file, to avoid an excessively large file
    // when there are many clumps
    constexpr int MaxPlottedClumps = 100;
}

//////////////////////////////////////////////////////////////////////

double ClumpySphericalSpatialGrid::meshEpsilonScale() const
{
    // leave at 1e-10 -- with 1e-11, inaccuracies start to occur.
    return 1e-10;
}

////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::setupSelfAfter()
{
    StructuredSphereSpatialGrid::setupSelfAfter();
    Log* log = find<Log>();

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
    _bvh = new SphericalClumpBVH;
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
            int s = index(i, j, k);

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
    // number of structured cells overlapped by each clump; only needed during setup, for logging
    vector<int> numOverlappingCellsPerClump(_numClumps, 0);
    log->info("Estimating structured-cell volumes overlapped by " + std::to_string(_numClumps) + " clumps (at least "
              + std::to_string(K) + " samples per clump)...");
    vector<std::pair<int, int>> hits;  // reused across clumps to avoid reallocating on every iteration
    for (int ci = 0; ci != _numClumps; ++ci)
    {
        const Clump& clump = _clumps[ci];
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
        numOverlappingCellsPerClump[ci] = numOverlappingCells;
        if (numOverlappingCells > 1) numStraddling++;

        for (const auto& hit : hits)
        {
            double fraction = static_cast<double>(hit.second) / totalSamples;
            _cellVolume[hit.first] -= clumpVolume * fraction;

            // if this cell is entirely inside the clump, its true remaining volume is exactly
            // zero regardless of what the Monte Carlo estimate above suggests; overriding it here
            // sidesteps any inaccuracy in that estimate for a cell much smaller than the clump,
            // and -- more importantly -- guarantees the exact zero that randomPositionInCell()
            // relies on to avoid an unbounded rejection loop
            double rmin, thetamin, phimin, rmax, thetamax, phimax;
            getCoords(hit.first, rmin, thetamin, phimin, rmax, thetamax, phimax);
            if (Quadrics::isSphericalCellInSphere(rmin, rmax, thetamin, thetamax, phimin, phimax, clump.center(),
                                                  clump.radius()))
                _cellVolume[hit.first] = 0.;
        }
    }

    // clamp away any negative volumes caused by Monte Carlo noise or near-coincident boundaries
    for (double& v : _cellVolume) v = max(v, 0.);

    // inform the user
    log->info("Summary:");
    log->info("  Clumps read from file:      " + std::to_string(_numClumps + numOutside + numOverlapping));
    log->info("  Rejected (outside domain):  " + std::to_string(numOutside));
    log->info("  Rejected (overlapping):     " + std::to_string(numOverlapping));
    log->info("  Remaining in spatial grid:  " + std::to_string(_numClumps));
    log->info("  Straddling a cell wall:     " + std::to_string(numStraddling));

    // report BVH leaf statistics, a diagnostic aid for tuning the tree build parameters
    // (NumBins and MaxLeafSize in SphericalClumpBVH.cpp); leaf diagonals are given relative to
    // the domain's outer radius to make them dimensionless and comparable across models
    ///auto leafStats = _bvh->leafStatistics();
    ///log->info("Statistics:");
    ///log->info("  Leaf depth:              min " + std::to_string(leafStats.minDepth) + ", max "
    ///          + std::to_string(leafStats.maxDepth) + ", avg " + StringUtils::toString(leafStats.avgDepth, 'f', 2));
    ///log->info("  Clumps per leaf:         min " + std::to_string(leafStats.minCount) + ", max "
    ///          + std::to_string(leafStats.maxCount) + ", avg " + StringUtils::toString(leafStats.avgCount, 'f', 2));
    ///log->info("  Leaf diagonal / rmax:    min " + StringUtils::toString(leafStats.minDiagonal / maxRadius(), 'f', 4)
    ///          + ", max " + StringUtils::toString(leafStats.maxDiagonal / maxRadius(), 'f', 4) + ", avg "
    ///          + StringUtils::toString(leafStats.avgDiagonal / maxRadius(), 'f', 4));

    // ---- epsilon used to progress a path, small relative to both the mesh and the smallest clump ----

    if (_numClumps)
    {
        double smallestClumpRadius = DBL_MAX;
        for (const Clump& clump : _clumps) smallestClumpRadius = min(smallestClumpRadius, clump.radius());
        _eps = min(_eps, 0.1 * smallestClumpRadius);
    }
}

//////////////////////////////////////////////////////////////////////

ClumpySphericalSpatialGrid::~ClumpySphericalSpatialGrid()
{
    delete _bvh;
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
    return StructuredSphereSpatialGrid::volume(m - _numClumps);
}

//////////////////////////////////////////////////////////////////////

double ClumpySphericalSpatialGrid::diagonal(int m) const
{
    if (m >= 0 && m < _numClumps) return 2. * _clumps[m].radius();
    return StructuredSphereSpatialGrid::diagonal(m - _numClumps);
}

//////////////////////////////////////////////////////////////////////

int ClumpySphericalSpatialGrid::cellIndex(Position bfr) const
{
    int c = _bvh->anyClumpContaining(bfr);
    if (c >= 0) return c;

    int s = StructuredSphereSpatialGrid::cellIndex(bfr);
    return s < 0 ? -1 : _numClumps + s;
}

//////////////////////////////////////////////////////////////////////

Position ClumpySphericalSpatialGrid::centralPositionInCell(int m) const
{
    if (m >= 0 && m < _numClumps) return _clumps[m].center();

    Position candidate = StructuredSphereSpatialGrid::centralPositionInCell(m - _numClumps);

    // if the nominal center happens to fall inside one of the clumps, fall back to a random position instead
    if (_bvh->anyClumpContaining(candidate) >= 0) return randomPositionInCell(m);
    return candidate;
}

//////////////////////////////////////////////////////////////////////

Position ClumpySphericalSpatialGrid::randomPositionInCell(int m) const
{
    if (m >= 0 && m < _numClumps) return random()->positionInSphere(_clumps[m].center(), _clumps[m].radius());

    int s = m - _numClumps;
    if (s < 0 || s >= _Ncells) return Position();

    // if this cell has no free volume left -- most commonly because it is entirely inside a
    // single clump, detected exactly during setup -- no rejection-sampled candidate could ever be
    // accepted, so fall back directly to the cell's nominal center instead of looping forever
    if (_cellVolume[s] <= 0.) return StructuredSphereSpatialGrid::centralPositionInCell(s);

    // safety net for the one configuration setup does not check for exactly: a cell fully covered
    // not by a single clump but by the union of several disjoint ones, leaving a nonzero but
    // vanishingly small (or, due to Monte Carlo noise, nonzero but spurious) free volume; cap the
    // number of rejection attempts and fall back to the same nominal center rather than hang
    constexpr int MaxAttempts = 1000;
    for (int attempt = 0; attempt != MaxAttempts; ++attempt)
    {
        Position candidate = StructuredSphereSpatialGrid::randomPositionInCell(s);
        if (_bvh->anyClumpContaining(candidate) < 0) return candidate;
    }
    return StructuredSphereSpatialGrid::centralPositionInCell(s);
}

//////////////////////////////////////////////////////////////////////

class ClumpySphericalSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const ClumpySphericalSpatialGrid* _grid{nullptr};
    double _eps{0.};
    int _clump{-1};              // index of the clump currently inside, or -1 for the structured grid
    int _i{-1}, _j{-1}, _k{-1};  // structured bin indices; meaningful only while not inside a clump

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
                if (_clump < 0 && !setCellIndices()) return abortPath();
                setState(State::Inside);

                // intentionally falls through to determine the first actual segment
            }

            // intentionally falls through
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
                        setSegment(_grid->_numClumps + _grid->index(icur, jcur, kcur), dsClump);
                        propagater(dsClump + _eps);
                        _clump = clumpHit;
                    }
                    else
                    {
                        setSegment(_grid->_numClumps + _grid->index(icur, jcur, kcur), ds);
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
        return false;
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
    // structured grid
    StructuredSphereSpatialGrid::write_xy(outfile);

    // spherical clumps crossing the xy plane
    for (const auto& hit : _bvh->clumpsCrossingPlane(2, 0.))
    {
        const Clump& c = _clumps[hit.first];
        outfile->writeCircle(c.center().x(), c.center().y(), hit.second);
    }
}

//////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    // structured grid: spheres and cones
    writeMeridionalStructure(outfile);

    // spherical clumps crossing the xz plane
    for (const auto& hit : _bvh->clumpsCrossingPlane(1, 0.))
    {
        const Clump& c = _clumps[hit.first];
        outfile->writeCircle(c.center().x(), c.center().z(), hit.second);
    }
}

//////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    // structured grid: spheres and cones -- identical to the xz view since the structured grid is axisymmetric
    writeMeridionalStructure(outfile);

    // spherical clumps crossing the yz plane
    for (const auto& hit : _bvh->clumpsCrossingPlane(0, 0.))
    {
        const Clump& c = _clumps[hit.first];
        outfile->writeCircle(c.center().y(), c.center().z(), hit.second);
    }
}

//////////////////////////////////////////////////////////////////////

void ClumpySphericalSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // structured grid
    StructuredSphereSpatialGrid::write_xyz(outfile);

    // a coarse wireframe sphere for each clump, up to a maximum to avoid an excessively large file
    int numPlotted = min(_numClumps, MaxPlottedClumps);
    for (int m = 0; m != numPlotted; ++m)
    {
        Position center = _clumps[m].center();
        outfile->writeSphere(center.x(), center.y(), center.z(), _clumps[m].radius());
    }
    if (_numClumps > MaxPlottedClumps)
    {
        find<Log>()->info("Plotted only the first " + std::to_string(MaxPlottedClumps) + " of "
                          + std::to_string(_numClumps) + " clumps in the 3D grid plot file");
    }
}

//////////////////////////////////////////////////////////////////////

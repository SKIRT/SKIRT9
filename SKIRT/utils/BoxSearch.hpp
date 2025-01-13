/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXSEARCH_HPP
#define BOXSEARCH_HPP

#include "Array.hpp"
#include "Box.hpp"
#include <functional>
#include <unordered_set>

//////////////////////////////////////////////////////////////////////

/** BoxSearch is a utility class for organizing spatial objects in a data structure that allows
    efficient retrieval of all objects that overlap a given point or ray. The class actually works
    with the bounding boxes of the objects being held, and leaves more detailed tests for
    containment or intersection to the client code.

    The spatial objects held by a BoxSearch instance are called entities. They are identified by a
    unique index \f$m\f$ ranging from 0 to \f$M-1\f$, where \f$M\f$ is the number of managed
    entities. All entities are handed to the BoxSearch instance in one go, so that they can be
    "bulk-loaded" into the search structure.

    The current implementation proceeds as follows. First, a regular Cartesian grid is contructed
    that partitions 3D space into \f$N_b^3\f$ blocks, where \f$N_b\f$ depends on the number of
    entities. Each block is then assigned a list of indices for all entities that possibly
    intersect with the block. In an attempt to balance the list lengths, the block separation
    points in each coordinate direction are chosen so that the entity bounding box centers are
    approximately evently distributed over the blocks in that direction. Locating the block
    containing a given query position then boils down to three binary searches (one in each
    direction). */
class BoxSearch
{
    // ------- Constructing and loading -------

public:
    /** The constructor creates a trivial BoxSearch instance holding no entities. Any queries
        will come up empty. */
    BoxSearch();

    /** This function loads the specified number of entities \f$M\f$ into the search structure,
        using the information returned by the provided callback functions. Any entities held
        prevously are removed and replaced by the new ones.

        The \em bounds callback function returns the bounding box of the entity with the given
        index. The \em intersects callback function returns true if the entity with the given index
        possibly intersects the specified box, and false otherwise. A relaxed intersection test
        such as testing against the bounding box is allowed but a stricter intersection test will
        result in more efficient retrieval queries.

        The callback functions are invoked one or more times for indices \f$m\f$ ranging from 0 to
        \f$M-1\f$, in arbitrary order. For given values of their argument(s), the callback
        functions must always return the same value. */
    void loadEntities(int numEntities, std::function<Box(int m)> bounds,
                      std::function<bool(int m, const Box& box)> intersects);

    // ------- Getting properties and statistics -------

public:
    /** This function returns the extent of the search domain, i.e. the union of all entity
        bounding boxes. */
    const Box& extent() const;

    /** This function returns the number of blocks in the search structure for each spatial
        dimension, or zero if no entities have been loaded. */
    int numBlocks() const;

    /** This function returns the smallest number of entity references in a search block. */
    int minEntitiesPerBlock() const;

    /** This function returns the largest number of entity references in a search block. */
    int maxEntitiesPerBlock() const;

    /** This function returns the average number of entity references in a search block. */
    double avgEntitiesPerBlock() const;

    // ------- Private helper functions -------

private:
    /** This function returns the linear index in the block list vector for the block with indices
        (i,j,k) in the three spatial directions. */
    int blockIndex(int i, int j, int k) const;

    /** This function returns the linear index in the block list vector for the block containing
        the given position. */
    int blockIndex(Vec bfr) const;

    // ------- Private helper classes -------

private:
    /** This typedef defines the generator return type of the entitiesFor function for a position.
        It represents an iterable sequence of integers. We use a private typedef to hide the actual
        type, which in the current implementation is simply a reference to a standard vector. */
    using EntityGeneratorForPosition = const vector<int>&;

    /** This typedef defines the generator return type of the entitiesFor function for a ray.
        It represents an iterable sequence of integers. We use a private typedef to hide the actual
        type, which in the current implementation is simply a copy of a standard unordered set. */
    using EntityGeneratorForRay = std::unordered_set<int>;

    // ------- Querying -------

public:
    /** This function returns an iterable sequence of indices \f$m\f$ of all entities that may
        overlap the specified position, sorted in increasing order. The sequence may be empty.

        The function guarantees that the sequence includes all entities whose bounding box overlaps
        the position. On the other hand, the sequence may contain entities whose bounding box does
        \em not overlap the position. */
    EntityGeneratorForPosition entitiesFor(Vec bfr) const;

    /** This function returns an iterable sequence of indices \f$m\f$ of all entities that may
        overlap the specified ray (starting point and direction), in arbitrary order. The sequence
        may be empty.

        The function guarantees that the sequence includes all entities whose bounding box overlaps
        the ray. On the other hand, the sequence may contain entities whose bounding box does \em
        not overlap the ray. */
    EntityGeneratorForRay entitiesFor(Vec bfr, Vec bfk) const;

    // ------- Data members -------

private:
    // search structure
    Box _extent;                   // the extent of the search domain, i.e. the union of all bounding boxes
    int _numBlocks{0};             // the number of grid blocks nb in each direction
    Array _xgrid, _ygrid, _zgrid;  // the nb+1 grid separation points for each spatial direction
    vector<vector<int>> _listv;    // the nb*nb*nb lists of indices for entities overlapping each block
    vector<int> _empty;            // vector that stays empty, used for returning empty generator

    // statistics
    int _minEntitiesPerBlock{0};
    int _maxEntitiesPerBlock{0};
    double _avgEntitiesPerBlock{0};
};

//////////////////////////////////////////////////////////////////////

#endif

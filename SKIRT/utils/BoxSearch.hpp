/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXSEARCH_HPP
#define BOXSEARCH_HPP

#include "Array.hpp"
#include "Box.hpp"
#include "Direction.hpp"
#include "EntityCollection.hpp"
#include "Position.hpp"

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
    entities with a floor of \f$N_b=20\f$ up to 1 million entities. Each block is then assigned a
    list of indices for all entities that possibly intersect with the block. In an attempt to
    balance the list lengths, the block separation points in each coordinate direction are chosen
    so that the entity bounding box centers are approximately evently distributed over the blocks
    in that direction. Locating the block containing a given query position then boils down to
    three binary searches (one in each direction).

    For performance reasons, the query functions use a template argument type instead of the
    appropriate std::function<> type declaration, which means they must be implemented inline in
    this header. */
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

    // ------- Querying -------

public:
    /** This function returns the index \f$m\f$ of the first entity containing the specified
        position, or -1 if there is no such entity. The callback function must return true if the
        entity with given index contains the position passed to the main function, and false
        otherwise. It is called zero or more times with the following guarantees: (1) all entities
        whose bounding box overlaps the position are included in the invocation sequence, and (2)
        the indices of consecutive invocations are sorted in increasing order.

        The first index for which the callback function returns true is returned to the caller of
        the main function. In case multiple entities overlap the position, because of the second
        guarantee, the entity with the smallest index is always returned.

        Note that the callback function may be called for entities whose bounding box does \em not
        overlap the position; in that case it should return false. */
    template<typename F> int firstEntity(Vec bfr, /* std::function<bool(int m)> */ F contains) const;

    /** This function accumulates the weights returned by the callback function for all entities
        that overlap the specified position, and returns the sum. If there are no such entities,
        the function returns zero.

        The callback function must return the weight in the entity with given index for the
        position passed to the main function, or zero if the entity does not contain the position.
        The callback function will be invoked with a sequence of indices in arbitrary order but
        without duplicates.

        Note that the callback function may be called for entities whose bounding box does \em not
        overlap the position; in that case it should return zero. */
    template<typename F> double accumulate(Vec bfr, /* std::function<double(int m)> */ F weight) const;

    /** This function replaces the contents of the specified entity collection by the set of
        entities that overlap the specified position with a nonzero weight. If there are no such
        entities, the collection will be empty.

        The callback function must return the weight in the entity with given index for the
        position passed to the main function, or zero if the entity does not contain the position.
        The callback function will be invoked with a sequence of indices in arbitrary order but
        without duplicates.

        Note that the callback function may be called for entities whose bounding box does \em not
        overlap the position; in that case it should return zero. */
    template<typename F>
    void getEntities(EntityCollection& entities, Vec bfr, /* std::function<double(int m) */ F weight) const;

    /** This function replaces the contents of the specified entity collection by the set of
        entities that overlap the specified ray (starting point and direction) with a nonzero
        weight. If there are no such entities, the collection will be empty.

        The callback function must return the weight in the entity with given index for the ray
        passed to the main function, or zero if the entity does not intersect the ray. The callback
        function will be invoked with a sequence of indices in arbitrary order and possibly with
        duplicates. For a given index, the callback function must always return the same value.
        Even if the callback function is invoked multiple times for the same index, the
        corresponding entity will be added to the entity collection just once.

        Note that the callback function may be called for entities whose bounding box does \em not
        overlap the position; in that case it should return zero. */
    template<typename F>
    void getEntities(EntityCollection& entities, Position bfr, Direction bfk,
                     /* std::function<double(int m)> */ F weight) const;

    // ------- Private helper functions -------

private:
    /** This function returns the linear index in the block list vector for the block with indices
        (i,j,k) in the three spatial directions. */
    int blockIndex(int i, int j, int k) const;

    /** This function returns the linear index in the block list vector for the block containing
        the given position. */
    int blockIndex(Vec bfr) const;

    /** If the specified ray (starting point and direction) intersects the search domain, this
        function returns true and stores the index ranges in three directions for the blocks
        overlapping the ray's bounding box in the arguments. If the ray does not intersect the
        domain, the function returns false and the contents of the argument indices is unspecified.
        */
    bool indexRangeForRay(Position bfr, Direction bfk, int& i1, int& i2, int& j1, int& j2, int& k1, int& k2) const;

    // ------- Data members -------

private:
    // search structure
    Box _extent;                   // the extent of the search domain, i.e. the union of all bounding boxes
    int _numBlocks{0};             // the number of grid blocks nb in each direction
    Array _xgrid, _ygrid, _zgrid;  // the nb+1 grid separation points for each spatial direction
    vector<vector<int>> _listv;    // the nb*nb*nb lists of indices for entities overlapping each block

    // statistics
    int _minEntitiesPerBlock{0};
    int _maxEntitiesPerBlock{0};
    double _avgEntitiesPerBlock{0};
};

//////////////////////////////////////////////////////////////////////

template<typename F> int BoxSearch::firstEntity(Vec bfr, F contains) const
{
    if (_numBlocks)
    {
        // loop over all entities overlapping the block containing the position
        for (int m : _listv[blockIndex(bfr)])
        {
            // return the first entity that actually contains the position
            if (contains(m)) return m;
        }
    }

    // there is no entity containing the position
    return -1;
}

//////////////////////////////////////////////////////////////////////

template<typename F> double BoxSearch::accumulate(Vec bfr, F weight) const
{
    double sum = 0.;

    if (_numBlocks)
    {
        // sum the weights for all entities overlapping the block containing the position
        for (int m : _listv[blockIndex(bfr)]) sum += weight(m);
    }
    return sum;
}

//////////////////////////////////////////////////////////////////////

template<typename F> void BoxSearch::getEntities(EntityCollection& entities, Vec bfr, F weight) const
{
    entities.clear();

    if (_numBlocks)
    {
        // add all entities overlapping that block to the collection with their respective weight
        for (int m : _listv[blockIndex(bfr)]) entities.add(m, weight(m));
    }
}

//////////////////////////////////////////////////////////////////////

template<typename F>
void BoxSearch::getEntities(EntityCollection& entities, Position bfr, Direction bfk, F weight) const
{
    entities.clear();

    if (_numBlocks)
    {
        // find the indices for first and last block, in each spatial direction,
        // overlapped by the bounding box of the path's intersection with the domain
        int i1, i2, j1, j2, k1, k2;
        if (indexRangeForRay(bfr, bfk, i1, i2, j1, j2, k1, k2))
        {
            // loop over all blocks in that 3D range
            for (int i = i1; i <= i2; i++)
                for (int j = j1; j <= j2; j++)
                    for (int k = k1; k <= k2; k++)
                    {
                        // if the path intersects the block
                        Box block(_xgrid[i], _ygrid[j], _zgrid[k], _xgrid[i + 1], _ygrid[j + 1], _zgrid[k + 1]);
                        double smin, smax;
                        if (block.intersects(bfr, bfk, smin, smax))
                        {
                            // add all entities overlapping that block to the collection with their respective weight
                            for (int m : _listv[blockIndex(i, j, k)])
                            {
                                entities.add(m, weight(m));
                            }
                        }
                    }
        }
    }
}

//////////////////////////////////////////////////////////////////////

#endif

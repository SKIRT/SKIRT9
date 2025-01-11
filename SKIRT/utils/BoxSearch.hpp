/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXSEARCH_HPP
#define BOXSEARCH_HPP

#include "Array.hpp"
#include "Box.hpp"
#include "Direction.hpp"
#include "Position.hpp"
#include <functional>
class EntityCollection;

//////////////////////////////////////////////////////////////////////

/** BoxSearch is a helper class for organizing spatial objects in a data structure that allows
    efficient retrieval of all objects that overlap a given point or ray. The class actually works
    with the bounding boxes of the objects being held, and leaves more detailed tests for
    containment or intersection to the client code.

    The spatial objects held by a BoxSearch instance are called entities. They are identified by
    a unique index \f$m\f$ ranging from 0 to \f$M-1\f$, where \f$M\f$ is the number of managed
    entities.

    To be completed. */
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
    int firstEntity(Vec bfr, std::function<bool(int m)> contains) const;

    /** This function accumulates the weights returned by the callback function for all entities
        that overlap the specified position, and returns the sum. If there are no such entities,
        the function returns zero.

        The callback function must return the weight in the entity with given index for the
        position passed to the main function, or zero if the entity does not contain the position.
        The callback function will be invoked with a sequence of indices in arbitrary order but
        without duplicates.

        Note that the callback function may be called for entities whose bounding box does \em not
        overlap the position; in that case it should return zero. */
    double accumulate(Vec bfr, std::function<double(int m)> weight) const;

    /** This function replaces the contents of the specified entity collection by the set of
        entities that overlap the specified position with a nonzero weight. If there are no such
        entities, the collection will be empty.

        The callback function must return the weight in the entity with given index for the
        position passed to the main function, or zero if the entity does not contain the position.
        The callback function will be invoked with a sequence of indices in arbitrary order but
        without duplicates.

        Note that the callback function may be called for entities whose bounding box does \em not
        overlap the position; in that case it should return zero. */
    void getEntities(EntityCollection& entities, Vec bfr, std::function<double(int m)> weight) const;

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
    void getEntities(EntityCollection& entities, Position bfr, Direction bfk,
                     std::function<double(int m)> weight) const;

    // ------- Private helper functions -------

private:
    /** This function returns the linear index in the block list vector for the block with indices
        (i,j,k) in the three spatial directions. */
    int blockIndex(int i, int j, int k) const;

    /** This function returns the linear index in the block list vector for the block containing
        the given position. */
    int blockIndex(Vec bfr) const;

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

#endif

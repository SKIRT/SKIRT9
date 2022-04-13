/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ENTITYCOLLECTION_HPP
#define ENTITYCOLLECTION_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** EntityCollection is a helper class used with the Snapshot class and its subclasses. As
    indicated there, some snapshot types use smoothed particles and some use spatial cells as their
    basic constituents. The generic term \em entity refers to either particles or cells. Assuming a
    snapshot with \f$N_\mathrm{ent}\f$ entities, any given entity is identfied by its zero-based
    index \f$0\le m \le N_\mathrm{ent}-1\f$.

    An EntityCollection object contains an unordered list of entity indices with corresponding
    weights. This can be used to represent the set of entities that overlap a given position or
    path in space with their corresponding relative weights. The relative weights
    would, depending on the circumstances, take into account the smoothing kernel of an imported
    particle and/or the length of the path segment overlapping the particle or cell.

    For performance reasons, the implementation allows the same entity to occur multiple times in
    the collection, although this is not recommended. Furthermore, an EntityCollection instance can
    be reused for storing a new collection without causing memory reallocations. */
class EntityCollection
{
public:
    // ------- Constructor -------

    /** The constructor creates an empty collection. */
    EntityCollection();

    // ------- Adding entities -------

    /** This function removes any entities from the collection, resulting in an empty collection. */
    void clear();

    /** This function adds an entity with index \f$m\f$ and weight \f$w\f$ to the collection. If
        \f$m<0\f$ or \f$w\le 0\f$, the function does nothing. */
    void add(int m, double w);

    /** This function removes any pre-existing entities from the collection and then adds a single
        entity with index \f$m\f$ and arbitrary weight 1 to the collection. If \f$m<0\f$, the
        function clears the collection but does not add an entity. */
    void addSingle(int m);

    // ------- Retrieving entities -------

    // basic data structure holding information about an entity
    struct Entity
    {
        int m;     // entity index
        double w;  // corresponding weight
    };

    /** This function returns (a read-only reference to) a list of the current entities in the
        collection, in arbitrary order. Each Entity object in the list holds two public data
        members: the index \em m of the entity being represented and the corresponding weight
        \f$w\f$. The entity index is guaranteed to be nonnegative. */
    const vector<Entity>& entities() const { return _entities; }

    /** This function returns true if the collection currently contains exactly one entity, and
        false otherwise. */
    bool isSingle() const { return _entities.size() == 1; }

    /** If the collection currently contains exactly one entity, i.e. isSingle() returns true, this
        function returns the index of that single entity. Otherwise, the behavior of this function
        is undefined. */
    bool singleIndex() const { return _entities[0].m; }

    // ------- Data members -------

private:
    vector<Entity> _entities;
};

//////////////////////////////////////////////////////////////////////

#endif

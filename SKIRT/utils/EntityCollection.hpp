/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ENTITYCOLLECTION_HPP
#define ENTITYCOLLECTION_HPP

#include "Basics.hpp"
#include <functional>
#include <unordered_map>
class Vec;

//////////////////////////////////////////////////////////////////////

/** EntityCollection is a helper class used with the Snapshot class and its subclasses. As
    indicated there, some snapshot types use smoothed particles and some use spatial cells as their
    basic constituents. The generic term \em entity refers to either particles or cells. Assuming a
    snapshot with \f$N_\mathrm{ent}\f$ entities, any given entity is identfied by its zero-based
    index \f$0\le m \le N_\mathrm{ent}-1\f$.

    An EntityCollection object contains an unordered collection of unique entity indices with
    corresponding weights. This can be used to represent the set of entities that overlap a given
    position or path in space with their corresponding relative weights. The relative weights
    would, depending on the circumstances, take into account the smoothing kernel of an imported
    particle and/or the length of the path segment overlapping the particle or cell.

    An EntityCollection instance can be reused for storing a new collection without forcing memory
    reallocations. */
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
        \f$m<0\f$, \f$w\le 0\f$, \f$w\f$ is NaN or infinity, or an entity with the same index is
        already in the collection, the function does nothing. */
    void add(int m, double w);

    /** This function removes any pre-existing entities from the collection and then adds a single
        entity with index \f$m\f$ and weight 1 to the collection. If \f$m<0\f$, the function clears
        the collection but does not add an entity. */
    void addSingle(int m);

    // ------- Retrieving entities -------

    /** Type declaration for the unordered map underlying the entity collection. */
    using Map = std::unordered_map<int, double>;

    /** This function returns read-only iterator to the first entity in the collection. Each entity
        consists of a key/value pair, where the key is the index \f$m\f$ of the entity being
        represented and the value is the corresponding weight \f$w\f$. The entity index is
        guaranteed to be nonnegative. */
    Map::const_iterator begin() const { return _entities.cbegin(); }

    /** This function returns read-only iterator just beyond the last entity in the collection. */
    Map::const_iterator end() const { return _entities.cend(); }

    // ------- Calculating results -------

    /** This function returns the accumulated weighted value of a given scalar field \f$f(m)\f$
        over all entities in the collection, \f[\sum_m f(m)\,w_m.\f] The argument specifies the
        scalar field \f$f(m)\f$; the call-back function should return the field value corresponding
        to a given entity index. */
    double accumulate(std::function<double(int m)> value);

    /** This function returns the nominator and denominator for the weighted average of a given
        scalar field \f$f(m)\f$ with given external weight \f$\omega(m)\f$ over all entities in the
        collection. The arguments respectively specify the scalar field \f$f(m)\f$ and the
        corresponding external weight \f$\omega(m)\f$. These functions should return the field
        value respectively the external weight corresponding to a given entity index.

        The function combines the external weights with the weights stored internally for each
        entity in the collection. Specifically, it calculates \f$\sum_m f(m)\,\omega(m)\,w_m\f$ and
        \f$\sum_m \omega(m)\,w_m\f$. */
    std::pair<double, double> average(std::function<double(int m)> value, std::function<double(int m)> weight);

    /** This function returns the nominator and denominator for the weighted average of a given
        vector field \f${\bf{f}}(m)\f$ with given external weight \f$\omega(m)\f$ over all entities
        in the collection. The arguments respectively specify the vector field \f${\bf{f}}(m)\f$
        and the corresponding external weight \f$\omega(m)\f$. These functions should return the
        field value respectively the external weight corresponding to a given entity index.

        The function combines the external weights with the weights stored internally for each
        entity in the collection. Specifically, it calculates \f$\sum_m
        {\bf{f}}(m)\,\omega(m)\,w_m\f$ and \f$\sum_m \omega(m)\,w_m\f$. */
    std::pair<Vec, double> average(std::function<Vec(int m)> value, std::function<double(int m)> weight);

    // ------- Data members -------

private:
    Map _entities;
};

//////////////////////////////////////////////////////////////////////

#endif

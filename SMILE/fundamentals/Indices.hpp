/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INDICES_HPP
#define INDICES_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The Indices class can be used to iterate over zero-based indices of type int
    in forward or backward order based on a run-time flag. For example:

    \code
    for (int index : Indices(4)) {  }         // loops over 0,1,2,3
    for (int index : Indices(4, true)) {  }   // loops over 3,2,1,0
    \endcode

    It is also possible to specify any container that implements size():

    \code
    vector xv(4);
    for (int index : Indices(xv)) {  }        // loops over 0,1,2,3
    for (int index : Indices(xv, true)) {  }  // loops over 3,2,1,0
    \endcode

    */
class Indices
{
private:
    /** Instances of this helper class are returned as index iterators. */
    class Iterator
    {
    public:
        Iterator(int value, int step) : _value(value), _step(step) {}
        bool operator!=(Iterator const& other) const { return _value != other._value; }
        int operator*() const { return _value; }
        Iterator& operator++()
        {
            _value += _step;
            return *this;
        }

    private:
        int _value;
        int _step;
    };

public:
    /** This constructor intializes an index range with the specified number of zero-based indices.
        The second argument indicates iteration direction: false (or missing) for forward
        iteration; true for backward iteration. */
    Indices(int range, bool reverse = false) : _range(range), _reverse(reverse) {}

    /** This constructor intializes an index range with a number of zero-based indices
        corresponding to the size of the specified container. The container class must implement
        the size() function. The second argument indicates iteration direction: false (or missing)
        for forward iteration; true for backward iteration. */
    template<class C> Indices(C cont, bool reverse = false) : _range(cont.size()), _reverse(reverse) {}

    /** This function returns the begin iterator corresponding to the iteration direction. */
    Iterator begin() const { return _reverse ? Iterator(_range - 1, -1) : Iterator(0, 1); }

    /** This function returns the end iterator corresponding to the iteration direction. */
    Iterator end() const { return _reverse ? Iterator(-1, -1) : Iterator(_range, 1); }

private:
    int _range;
    bool _reverse;
};

////////////////////////////////////////////////////////////////////

#endif

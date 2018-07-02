/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RANGE_HPP
#define RANGE_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** Range is a low-level class for working with ranges: each instance represents a range with a
    minimum and maximum value. The class is fully implemented inline (in this header file). Most
    compilers optimize away all overhead so that using this class is just as efficient as directly
    writing the code in terms of the range components. */
class Range
{
protected:
    double _min, _max;

public:
    /** This is the default constructor; the minimum and maximum values are initialized to zero. */
    Range() : _min(0), _max(0) { }

    /** This constructor initializes the range to the minimum and maximum values provided as arguments. */
    Range(double min, double max) : _min(min), _max(max) { }

    /** This function sets the range to the minimum and maximum values provided as arguments. */
    void set(double min, double max) { _min = min; _max = max; }

    /** This function returns the minimum value of the range. */
    double min() const { return _min; }

    /** This function returns the maximum value of the range. */
    double max() const { return _max; }

    /** This function returns the middle value of the range. */
    double mid() const { return 0.5 * (_min + _max); }

    /** This function returns the width of the range, i.e. the maximum value minus the minimum value. */
    double width() const { return _max - _min; }

    /** This function returns true if the given value is inside the range, and false otherwise. */
    bool contains(double x) const { return x >= _min && x <= _max; }

    /** This function returns true if the given value is inside the range, with the given fuzzyness
        factor, and false otherwise. */
    bool containsFuzzy(double x, double eps=1e-14) const { return x >= _min*(1-eps) && x <= _max*(1+eps); }

    /** This function returns true if the range is empty (its minimum is larger than or equal to
        its maximum), and false otherwise. */
    bool empty() const { return _min >= _max; }

    /** This function updates the range so that it represents the intersection of the original
        range with the other range given as an argument. If the two ranges do not overlap, the
        resulting range will have a minimum larger than or equal to its maximum. */
    void intersect(const Range& range)
    {
        if (_min < range._min) _min = range._min;
        if (_max > range._max) _max = range._max;
    }
};

//////////////////////////////////////////////////////////////////////

#endif

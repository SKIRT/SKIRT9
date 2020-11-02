/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHORTARRAY_HPP
#define SHORTARRAY_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** A ShortArray instance holds a an array of double values that is usually shorter than the
    compile-time constant ShortArray::N but can be larger if needed. As long as the array is not
    longer than N, no heap allocations occur. If the array does become longer than N, the
    appropriate space is allocated on the heap. The number of items held by the array can be
    adjusted, but only at the cost of erasing all previously stored values: the resize operation
    sets all values to zero, just as if the array was freshly constructed. ShortArray instances
    cannot be copied or moved.

    In the current implementation, the value of N is 8. */
class ShortArray
{
public:
    const static int N = 8;

    // ================== Constructing and assigning ==================

    /** The default constructor creates an empty array. */
    ShortArray() {}

    /** This constructor creates an empty array with the specified number of values, all
        initialized to zero. */
    explicit ShortArray(size_t n) : _n{n}
    {
        if (n > N)
        {
            _p = new double[n];
            _c = n;
        }
        clear();
    }

    /** The copy constructor is deleted. */
    ShortArray(const ShortArray&) = delete;

    /** The assignment operator is deleted. */
    ShortArray& operator=(const ShortArray&) = delete;

    /** The destructor releases the memory buffer if it was allocated on the heap. */
    ~ShortArray()
    {
        if (_c) delete[] _p;
    }

    // ================== Sizing and clearing ==================

    /** This function resizes the array so that it holds the specified number of values. The values
        are initialized to zero; the previous contents is lost. All iterators are invalidated. */
    void resize(size_t n)
    {
        _n = n;

        // the heap-allocated buffer, if it exists, always has a size larger than N,
        // so if the new size is smaller than N we're always safe regardless of the value of _c
        if (n > N && n > _c)
        {
            if (_c) delete[] _p;
            _p = new double[n];
            _c = n;
        }
        clear();
    }

    /** This function clears all values held by the array to zero. */
    void clear() { std::fill(_p, _p + _n, 0.); }

    /** This function returns the number of values held by the array. */
    size_t size() const { return _n; }

    // ================== Element access ==================

    /** This function returns a read-only reference to the specified array element. */
    double operator[](size_t i) const { return _p[i]; }

    /** This function returns a writable reference to the specified array element. */
    double& operator[](size_t i) { return _p[i]; }

private:
    // ================== Data members ==================

    double* _p{_a};  // points to the inside buffer if _c==0; to the heap-allocated buffer if _c>0
    size_t _n{0};    // the number of values held by the array
    size_t _c{0};    // the capacity of the heap-allocated buffer, or zero if there is none
    double _a[N];    // the inside buffer
};

////////////////////////////////////////////////////////////////////

#endif

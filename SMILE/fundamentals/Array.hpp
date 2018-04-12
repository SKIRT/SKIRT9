/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ARRAY_HPP
#define ARRAY_HPP

#include "Basics.hpp"
#include <valarray>

////////////////////////////////////////////////////////////////////

/** \class Array
An Array instance holds an array of double values, and allows easily and efficiently performing
mathematical operations on the corresponding values in multiple arrays. The number of items
held by the array can be adjusted, but only at the cost of erasing all previously stored
values: the resize operation sets all values to zero, just as if the array was freshly
constructed.

The copy/move constructors and assignment operators automatically resize the receiving array to
the size of the array being copied/moved. For optimal performance, the computed assignment
operators and the unary and binary arithmetic operators assume that all arrays involved in the
calculation have the same size. If not, the resulting behavior is undefined.

Below is a synopsis of the facilities offered by the Array class. In the current
implementation, Array is a type alias for \c std::valarray<double>. The \c std::valarray
template offers a number of facilities beyond those listed in the synopsis. To allow future
implementations to divert from std::valarray if needed, those extra facilities should not be
used. Specifically, do not use other constructors than those listed in the synopsis, and do not
use any of the helper classes (e.g. \c slice, \c gslice, \c mask_array, \c indirect_array).

\verbatim
class Array
{
public:
    // constructors
    Array();                            // constructs empty array
    explicit Array(size_t n);           // sets all values to zero
    Array(const Array& v);
    Array(Array&& v);
    ~Array();

    // assignment
    Array& operator=(const Array& v);
    Array& operator=(Array&& v);
    Array& operator=(double x);

    // swapping
    void swap(Array& v);

    // clearing and sizing
    void resize(size_t n);              // sets all values to zero after resize
    size_t size() const;

    // element access
    double operator[](size_t i) const;
    double& operator[](size_t i);

    // computed assignment - assume same size (undefined behavior if not)
    Array& operator+= (const Array& v);
    Array& operator-= (const Array& v);
    Array& operator*= (const Array& v);
    Array& operator/= (const Array& v);
    Array& operator+= (double x);
    Array& operator-= (double x);
    Array& operator*= (double x);
    Array& operator/= (double x);

    // unary operators
    Array operator+() const;
    Array operator-() const;

    // basic algorithms
    double sum() const;                 // undefined behavior for empty array
    double min() const;                 // undefined behavior for empty array
    double max() const;                 // undefined behavior for empty array

    // applying functions
    Array apply(double f(double)) const;
};

// swapping
void swap(Array& x, Array& y);

// binary operators - assume same size (undefined behavior if not)
Array operator+ (const Array& x, const Array& y);
Array operator- (const Array& x, const Array& y);
Array operator* (const Array& x, const Array& y);
Array operator/ (const Array& x, const Array& y);
Array operator+ (const Array& x, double y);
Array operator- (const Array& x, double y);
Array operator* (const Array& x, double y);
Array operator/ (const Array& x, double y);
Array operator+ (double x, const Array& y);
Array operator- (double x, const Array& y);
Array operator* (double x, const Array& y);
Array operator/ (double x, const Array& y);

// applying functions
Array abs(const Array& x);
Array acos(const Array& x);
Array asin(const Array& x);
Array atan(const Array& x);
Array atan2(const Array& y, const Array& x);
Array atan2(const Array& y, double x);
Array atan2(double y, const Array& x);
Array cos(const Array& x);
Array cosh(const Array& x);
Array exp(const Array& x);
Array log(const Array& x);
Array log10(const Array& x);
Array pow(const Array& x, const Array& y);
Array pow(const Array& x, double y);
Array pow(double x, const Array& y);
Array sin(const Array& x);
Array sinh(const Array& x);
Array sqrt(const Array& x);
Array tan(const Array& x);
Array tanh(const Array& x);

// iterators
double* begin(Array& v);
const double* begin(const Array& v);
double* end(Array& v);
const double* end(const Array& v);
\endverbatim
*/
using Array = std::valarray<double>;

////////////////////////////////////////////////////////////////////

#endif

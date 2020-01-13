/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VEC_HPP
#define VEC_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** Vec is a low-level class for working with three-dimensional vectors: each instance represents a
    vector with three cartesian components called \em x, \em y, and \em z. The class defines
    operators for adding two vectors and for multiplying a vector by a scalar. It offers functions
    to retrieve the vector components, to get the norm or the squared norm of a vector, and for
    calculating the dot product and cross product of two vectors. The Vec class is fully
    implemented inline (in this header file). Most compilers optimize away all overhead so that
    using this class is just as efficient as directly writing the code in terms of the vector
    components. */
class Vec
{
protected:
    /** These data members represent the cartesian vector components */
    double _x, _y, _z;

public:
    /** This is the default constructor; all vector components are initialized to zero. */
    inline Vec() : _x(0), _y(0), _z(0) {}

    /** This constructor initializes the vector components to the values provided as arguments. */
    inline Vec(double x, double y, double z) : _x(x), _y(y), _z(z) {}

    /** This function sets the vector components to the values provided as arguments. */
    inline void set(double x, double y, double z)
    {
        _x = x;
        _y = y;
        _z = z;
    }

    /** This function returns true if all components of the vector are trivially zero, and false
        otherwise. */
    inline bool isNull() const { return _x == 0. && _y == 0. && _z == 0.; }

    /** This function returns the \em x component of the vector. */
    inline double x() const { return _x; }

    /** This function returns the \em y component of the vector. */
    inline double y() const { return _y; }

    /** This function returns the \em z component of the vector. */
    inline double z() const { return _z; }

    /** This function returns the norm of the vector. */
    inline double norm() const { return sqrt(_x * _x + _y * _y + _z * _z); }

    /** This function returns the squared norm of the vector. */
    inline double norm2() const { return _x * _x + _y * _y + _z * _z; }

    /** This static function returns the dot product (inner product) of two vectors. */
    inline static double dot(Vec a, Vec b) { return a._x * b._x + a._y * b._y + a._z * b._z; }

    /** This static function returns the vector product (outer product) of two vectors. */
    inline static Vec cross(Vec a, Vec b)
    {
        return Vec(a._y * b._z - a._z * b._y, a._z * b._x - a._x * b._z, a._x * b._y - a._y * b._x);
    }

    /** This operator adds another vector to this one. */
    inline Vec& operator+=(Vec v)
    {
        _x += v._x;
        _y += v._y;
        _z += v._z;
        return *this;
    }

    /** This operator subtracts another vector from this one. */
    inline Vec& operator-=(Vec v)
    {
        _x -= v._x;
        _y -= v._y;
        _z -= v._z;
        return *this;
    }

    /** This operator multiplies this vector by a scalar. */
    inline Vec& operator*=(double s)
    {
        _x *= s;
        _y *= s;
        _z *= s;
        return *this;
    }

    /** This operator divides this vector by a scalar. */
    inline Vec& operator/=(double s)
    {
        _x /= s;
        _y /= s;
        _z /= s;
        return *this;
    }
};

//////////////////////////////////////////////////////////////////////

/** This free operator adds two vectors. */
inline Vec operator+(Vec a, Vec b)
{
    return a += b;
}

/** This free operator subtracts two vectors. */
inline Vec operator-(Vec a, Vec b)
{
    return a -= b;
}

/** This free operator multiplies a vector by a scalar. */
inline Vec operator*(Vec a, double s)
{
    return a *= s;
}

/** This free operator multiplies a scalar by a vector. */
inline Vec operator*(double s, Vec b)
{
    return b *= s;
}

/** This free operator divides a vector by a scalar. */
inline Vec operator/(Vec a, double s)
{
    return a /= s;
}

/** This free operator defines a simple ordering on vectors using the x coordinate. */
inline bool operator<(Vec a, Vec b)
{
    return a.x() < b.x();
}

//////////////////////////////////////////////////////////////////////

#endif

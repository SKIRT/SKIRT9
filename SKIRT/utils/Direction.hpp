/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DIRECTION_HPP
#define DIRECTION_HPP

#include "Vec.hpp"

//////////////////////////////////////////////////////////////////////

/** An object of the Direction class represents a direction in three-dimensional space. In
    principle only two coordinates are necessary to fully characterize a direction; a possible
    choice is spherical coordinates on the unit sphere. However in this class, a direction is
    internally represented by means of its three cartesian coordinates \f$(k_x,k_y,k_z)\f$. These
    coordinates must satisfy the normalization relation \f[ k_x^2+k_y^2+k_z^2=1. \f]

    To assist client code using the Direction class with ensuring this normalization, the
    constructors and the set() functions offer the Boolean \em normalize argument. If set to true,
    the constructor or function will automatically normalize the components passed by the caller.
    If set to false, the components are copied as given without any further verification. This can
    be used to avoid redundant calculations in case the caller has already ensured normalization.

    An important exception to the normalization invariant occurs when all three direction
    components are zero, resulting in an invalid "null" direction. This can happen, for example,
    when constructing a direction from user input or from the cross product of two parallel
    vectors. In situations where a null direction may occur, client code should test for this using
    the Vec::isNull() function.

    The Direction class is publicly based on the Vec class, and the direction coordinates
    \f$(k_x,k_y,k_z)\f$ are stored in the Vec components \em x, \em y and \em z. Consequently the
    functions and operators defined for Vec are automatically available for Direction, except for
    functions that would break the normalization (e.g. the compound-assignment operators). */
class Direction : public Vec
{
public:
    /** The default constructor for the class Direction creates a direction parallel to the
        \f$z\f$-axis of the reference frame. The cartesian coordinates of this direction are
        \f$(k_x,k_y,k_z)=(0,0,1)\f$. */
    inline Direction() : Vec(0, 0, 1) {}

    /** Constructor for a direction with given cartesian coordinates \f$k_x\f$, \f$k_y\f$ and
        \f$k_z\f$. If \em normalize is true, the constructor scales these coordinates so that they
        conform to the normalization relation \f$k_x^2+k_y^2+k_z^2=1\f$; if this is not possible
        because the coordinates are zero, the direction is set to the invalid null vector. If \em
        normalize is false, the components are copied as given without any further verification.
        This can be used to avoid redundant calculations in case the caller has already ensured
        normalization. */
    inline Direction(double kx, double ky, double kz, bool normalize) : Vec(kx, ky, kz)
    {
        if (normalize) this->normalize();
    }

    /** Constructor for a direction with given cartesian coordinates \f${\bf{k}}_x\f$,
        \f${\bf{k}}_y\f$ and \f${\bf{k}}_z\f$. If \em normalize is true, the constructor scales
        these coordinates so that they conform to the normalization relation
        \f$k_x^2+k_y^2+k_z^2=1\f$; if this is not possible because the coordinates are zero, the
        direction is set to the invalid null vector. If \em normalize is false, the components are
        copied as given without any further verification. This can be used to avoid redundant
        calculations in case the caller has already ensured normalization. */
    inline Direction(Vec k, bool normalize) : Vec(k)
    {
        if (normalize) this->normalize();
    }

    /** Constructor for a direction with a given azimuth and polar angle. The cartesian coordinates
        are calculated as \f[ \begin{split} k_x &= \sin\theta\,\cos\varphi, \\ k_y &=
        \sin\theta\,\sin\varphi, \\ k_z &= \cos\theta. \end{split} \f] It is explicitly checked
        whether \f$\theta\f$ lies within the interval \f$[0,\pi]\f$. If not so, this does not
        correspond to a viable direction, and a FatalError is thrown. */
    Direction(double theta, double phi);

    /** Disable the setter inherited from Vec so that callers are forced to specify a normalization
        option. */
    inline void set(double, double, double) = delete;

    /** This function sets the direction to the given cartesian coordinates \f${\bf{k}}_x\f$,
        \f${\bf{k}}_y\f$ and \f${\bf{k}}_z\f$. If \em normalize is true, the constructor scales
        these coordinates so that they conform to the normalization relation
        \f$k_x^2+k_y^2+k_z^2=1\f$; if this is not possible because the coordinates are zero, the
        direction is set to the invalid null vector. If \em normalize is false, the components are
        copied as given without any further verification. This can be used to avoid redundant
        calculations in case the caller has already ensured normalization. */
    inline void set(double kx, double ky, double kz, bool normalize)
    {
        _x = kx;
        _y = ky;
        _z = kz;
        if (normalize) this->normalize();
    }

    /** This function sets the direction to the given cartesian coordinates \f${\bf{k}}_x\f$,
        \f${\bf{k}}_y\f$ and \f${\bf{k}}_z\f$. If \em normalize is true, the constructor scales
        these coordinates so that they conform to the normalization relation
        \f$k_x^2+k_y^2+k_z^2=1\f$; if this is not possible because the coordinates are zero, the
        direction is set to the invalid null vector. If \em normalize is false, the components are
        copied as given without any further verification. This can be used to avoid redundant
        calculations in case the caller has already ensured normalization. */
    inline void set(Vec k, bool normalize)
    {
        _x = k.x();
        _y = k.y();
        _z = k.z();
        if (normalize) this->normalize();
    }

    /** This function returns the direction coordinate \f$k_x\f$. It is equivalent to the function
        x() inherited from the Vec class. */
    inline double kx() const { return _x; }

    /** This function returns the direction coordinate \f$k_y\f$. It is equivalent to the function
        y() inherited from the Vec class. */
    inline double ky() const { return _y; }

    /** This function returns the direction coordinate \f$k_z\f$. It is equivalent to the function
        z() inherited from the Vec class. */
    inline double kz() const { return _z; }

    /** This function returns the three cartesian coordinates \f$(k_x,k_y,k_z)\f$ of the direction.
        Note that these cartesian coordinates are also directly accessible through the functions
        kx(), ky() and kz(), and the equivalents x(), y() and z() inherited from the Vec class. */
    void cartesian(double& kx, double& ky, double& kz) const;

    /** This function determines the two spherical coordinates \f$(\theta,\varphi)\f$ of the
        direction. They are calculated from the internally stored cartesian coordinates by means of
        the formulae \f[ \begin{split} \theta &= \arccos\left(k_z\right), \\ \varphi &=
        \arctan\left(\frac{k_y}{k_x}\right). \end{split} \f] */
    void spherical(double& theta, double& phi) const;

    /** This operator returns the opposite direction (negating each coordinate). */
    Direction operator-() const { return Direction(-_x, -_y, -_z, false); }

    /** Disable this compound assignment operator inherited from Vec because it breaks the
        normalization. */
    inline Vec& operator+=(Vec v) = delete;

    /** Disable this compound assignment operator inherited from Vec because it breaks the
        normalization. */
    inline Vec& operator-=(Vec v) = delete;

    /** Disable this compound assignment operator inherited from Vec because it breaks the
        normalization. */
    inline Vec& operator*=(Vec v) = delete;

    /** Disable this compound assignment operator inherited from Vec because it breaks the
        normalization. */
    inline Vec& operator/=(Vec v) = delete;

private:
    /** This private function normalizes the coordinates already stored in this instance by
        dividing them by the value of the Vec::norm() function. As an exception, if the norm()
        function returns zero, all coordinates are set to zero, resulting in an invalid direction.
        The function is called from constructors and setters for this class in case the client code
        sets \em normalize to true. */
    void normalize();
};

//////////////////////////////////////////////////////////////////////

#endif

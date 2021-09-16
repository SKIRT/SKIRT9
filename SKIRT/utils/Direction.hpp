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
    coordinates must satisfy the normalization relation \f[ k_x^2+k_y^2+k_z^2=1. \f] Direction is
    publicly based on the Vec class, and the direction coordinates \f$(k_x,k_y,k_z)\f$ are stored
    in the Vec components \em x, \em y and \em z. Consequently all functions and operators
    defined for Vec are automatically available for Position as well. However this also means that
    the normalization relation can be broken by careless access to functions such as set()
    inherited from the Vec class. */
class Direction : public Vec
{
public:
    /** The default constructor for the class Direction creates a direction parallel to the
        \f$z\f$-axis of the reference frame. The cartesian coordinates of this direction are
        \f$(k_x,k_y,k_z)=(0,0,1)\f$. */
    inline Direction() : Vec(0, 0, 1) {}

    /** Constructor for a direction with given cartesian coordinates \f$k_x\f$, \f$k_y\f$ and
        \f$k_z\f$. In order to avoid abundant calculations, it is not checked whether the three
        cartesian coordinates satisfy the condition \f[ k_x^2+k_y^2+k_z^2=1. \f] It is the user's
        responsibility to provide valid coordinates. */
    inline Direction(double kx, double ky, double kz) : Vec(kx, ky, kz) {}

    /** Constructor for a direction with given cartesian coordinates \f${\bf{k}}_x\f$,
        \f${\bf{k}}_y\f$ and \f${\bf{k}}_z\f$. It is declared <tt>explicit</tt> to avoid implicit
        type conversions. In order to avoid abundant calculations, it is not checked whether the
        three cartesian coordinates satisfy the condition \f[ k_x^2+k_y^2+k_z^2=1. \f] It is the
        user's responsibility to provide valid coordinates. */
    explicit inline Direction(Vec k) : Vec(k) {}

    /** Constructor for a direction with a given azimuth and polar angle. The cartesian coordinates
        are calculated as \f[ \begin{split} k_x &= \sin\theta\,\cos\varphi, \\ k_y &=
        \sin\theta\,\sin\varphi, \\ k_z &= \cos\theta. \end{split} \f] It is explicitly checked
        whether \f$\theta\f$ lies within the interval \f$[0,\pi]\f$. If not so, this does not
        correspond to a viable direction, and a FatalError is thrown. */
    Direction(double theta, double phi);

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

    /** This operator returns the opposite direction (negating each component). */
    Direction operator-() const { return Direction(-_x, -_y, -_z); }
};

//////////////////////////////////////////////////////////////////////

#endif

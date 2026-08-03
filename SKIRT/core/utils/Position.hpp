/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef POSITION_HPP
#define POSITION_HPP

#include "Direction.hpp"
#include "Vec.hpp"

//////////////////////////////////////////////////////////////////////

/** An object of the Position class is used to define a position in three-dimensional space. To
    fully specify a position, three independent coordinates are necessary. Various options are
    possible. However in this class, a position is internally represented by means of its three
    cartesian coordinates \f$(x,y,z)\f$. Position is in fact publicly based on the Vec class, so
    that all functions and operators defined for Vec are automatically available for Position as
    well. */
class Position : public Vec
{
public:
    /** The default constructor for the class Position creates a point in the origin of the
        reference system, with cartesian coordinates \f$(x,y,z)=(0,0,0)\f$. */
    inline Position() : Vec() {}

    /** Constructor for a position with given cartesian coordinates \f$x\f$, \f$y\f$ and \f$z\f$.
        */
    inline Position(double x, double y, double z) : Vec(x, y, z) {}

    /** Constructor for a position with given cartesian coordinates \f${\bf{r}}_x\f$,
        \f${\bf{r}}_y\f$ and \f${\bf{r}}_z\f$. It is declared <tt>explicit</tt> to avoid implicit
        type conversions. */
    explicit inline Position(Vec r) : Vec(r) {}

    /** This enum has a constant for each of the supported coordinate systems. */
    enum class CoordinateSystem { CARTESIAN, CYLINDRICAL, SPHERICAL };

    /** Constructor for a position with three coordinates \f$(u,v,w)\f$ given in a coordinate
        system defined by the last argument. If the coordinate system is cartesian, the coordinates
        are just copied. If it is cylindrical, it is assumed that \f$(u,v,w) = (R,\phi,z)\f$ and
        the cartesian coordinates are calculated as \f[ \begin{split} x &=
        r\,\sin\theta\,\cos\varphi, \\ y &= r\,\sin\theta\,\sin\varphi, \\ z &= r\,\cos\theta,
        \end{split} \f] If spherical coordinates are used, we assume that \f$(u,v,w) =
        (r,\theta,\phi)\f$ and the cartesian coordinates are calculated as \f[ \begin{split} x &=
        R\,\cos\varphi, \\ y &= R\,\sin\varphi, \\ z &= z. \end{split} \f] */
    Position(double u, double v, double w, CoordinateSystem coordtype);

    /** Constructor for a position, starting from a radius \f$r\f$ and a direction \f${\bf k}\f$.
        With \f$(k_x,k_y,k_z)\f$ the cartesian coordinates of \f${\bf k}\f$, the construction of
        the cartesian coordinates for the position is straightforward, \f[ \begin{split} x &=
        r\,k_x, \\ y &= r\,k_y, \\ z &= r\,k_z. \end{split} \f] */
    Position(double r, Direction bfk);

    /** Constructor for a position, starting from only a direction \f${\bf k}\f$. This constructor
        hence generates a position on the unit sphere, and can be used to convert directions into
        positions. It is declared <tt>explicit</tt> to avoid implicit type conversions. */
    explicit Position(Direction bfk);

    /** This function returns the radius \f$r\f$ of the position, i.e. the distance between the
        position and the centre of the reference system, given by \f[ r=\sqrt{x^2+y^2+z^2}. \f]
        This function returns the result of the norm() function inherited from the Vec class. */
    double radius() const;

    /** This function returns the radius of the projection of the position on the XY-plane, which
        is the radial coordinate in a cylindrical coordinate system. The formula simply reads \f[ R
        = \sqrt{x^2+y^2}. \f] Useful in systems with an axial symmetry, where only \f$R\f$ and
        \f$z\f$ are necessary to characterize the physical conditions at a given position. */
    double cylRadius() const;

    /** This function returns the height of a position above the XY-plane, i.e. the
        \f$z\f$-coordinate. Useful in systems with an axial symmetry, where only \f$R\f$ and
        \f$z\f$ are necessary to characterize the physical conditions at a given position. */
    double height() const;

    /** This function returns the three cartesian coordinates \f$(x,y,z)\f$ of the position. Note
        that these cartesian coordinates are also directly accessible through the functions
        x(), y() and z() inherited from the Vec class. */
    void cartesian(double& x, double& y, double& z) const;

    /** This function determines the three spherical coordinates \f$(r,\theta,\phi)\f$ of the
        position. The connection between these coordinates and the internally stored cartesian
        coordinates is \f[ \begin{split} r &= \sqrt{x^2+y^2+z^2}, \\ \theta &=
        \arccos\left(\frac{z}{r}\right), \\ \varphi &= \arctan\left(\frac{y}{x}\right). \end{split}
        \f] */
    void spherical(double& r, double& theta, double& phi) const;

    /** This function determines the three cylindrical coordinates \f$(R,\varphi,z)\f$ of the
        position. The connection between these coordinates and the internally stored cartesian
        coordinates is \f[ \begin{split} R &= \sqrt{x^2+y^2}, \\ \varphi &=
        \arctan\left(\frac{y}{x}\right), \\ z &= z. \end{split} \f] */
    void cylindrical(double& R, double& phi, double& z) const;
};

//////////////////////////////////////////////////////////////////////

#endif

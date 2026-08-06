/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HOMOGENEOUSTRANSFORM_HPP
#define HOMOGENEOUSTRANSFORM_HPP

#include "Vec.hpp"

//////////////////////////////////////////////////////////////////////

/** Consider homogeneous coordinates \f$(x,y,z,w)\f$ with arbitrary scale factor
    \f$w\f$. An arbitary point \f$\mathbf{P}=(P_{x},P_{y},P_{z})\f$ can be represented
    in homogeneous coordinates as \f$\mathcal{P}=(P_{x},P_{y},P_{z},1)\f$.
    An arbitrary direction \f$\mathbf{D}=(D_{x},D_{y},D_{z})\f$ is represented
    as a point at infinity \f$\mathcal{D}=(D_{x},D_{y},D_{z},0)\f$.

    Vice versa, given the homogeneous coordinates \f$\mathcal{Q}=(Q_{x},Q_{y},Q_{z},Q_{w})\f$,
    the corresponding regular coordinates for the point can be retrieved
    as \f$\mathbf{Q}=(Q_{x}/Q_{w},Q_{y}/Q_{w},Q_{z}/Q_{w})\f$ provided \f$Q_{w}\neq0\f$.

    A homogeneous coordinate transformation can be written as \f$\mathcal{P}'=\mathcal{P}\,\mathcal{T}\f$
    where \f$\mathcal{T}\f$ is a \f$4\times4\f$ transformation matrix which
    can represent a scaling, rotation, translation, perspective transformation
    or any combination thereof; or more explicitly
    \f[
    \left[\begin{array}{cccc}
    x' & y' & z' & w'\end{array}\right]=\left[\begin{array}{cccc}
    x & y & z & w\end{array}\right]\left[\begin{array}{cccc}
    t_{11} & t_{12} & t_{13} & t_{14}\\
    t_{21} & t_{22} & t_{23} & t_{24}\\
    t_{31} & t_{32} & t_{33} & t_{34}\\
    t_{41} & t_{42} & t_{43} & t_{44}
    \end{array}\right]
    \f]

    A HomogeneousTransform instance represents the \f$4\times4\f$ matrix specifying a
    homogeneous coordinate transform. The constructor creates an identity transform, and
    various functions serve to concatenate specific transforms to the current one.
    Finally, the transform() function actually computes the transformation for given
    homogeneous coordinates. */
class HomogeneousTransform
{
public:
    /** The default constructor; it creates an identity transform. */
    HomogeneousTransform();

    /** This function concatenates to the current transform a translation over the specified offset
        along each coordinate axis. */
    void translate(double x, double y, double z);

    /** This function concatenates to the current transform a scaling by the specified factor
        along each coordinate axis. */
    void scale(double x, double y, double z);

    /** This function concatenates to the current transform a rotation about the x axis over an angle
        \f$\alpha\f$ specified through the values \f$\cos\alpha\f$ and \f$\sin\alpha\f$.
        It is the caller's responsibility to ensure that \f$\cos^2\alpha+\sin^2\alpha=1\f$
        holds for the values passed to this function.*/
    void rotateX(double cos, double sin);

    /** This function concatenates to the current transform a rotation about the y axis over an angle
        \f$\alpha\f$ specified through the values \f$\cos\alpha\f$ and \f$\sin\alpha\f$.
        It is the caller's responsibility to ensure that \f$\cos^2\alpha+\sin^2\alpha=1\f$
        holds for the values passed to this function.*/
    void rotateY(double cos, double sin);

    /** This function concatenates to the current transform a rotation about the z axis over an angle
        \f$\alpha\f$ specified through the values \f$\cos\alpha\f$ and \f$\sin\alpha\f$.
        It is the caller's responsibility to ensure that \f$\cos^2\alpha+\sin^2\alpha=1\f$
        holds for the values passed to this function.*/
    void rotateZ(double cos, double sin);

    /** This function concatenates to the current transform a perspective along the z axis
        with the specified focal length. */
    void perspectiveZ(double f);

    /** This function concatenates the specified tranform to the current transform. */
    void concatenate(const HomogeneousTransform& transform);

    /** This function applies the current transform to the specified homogeneous coordinates
        and stores the result in homogeneous coordinates in the output variables. */
    void transform(double x, double y, double z, double w, double& outx, double& outy, double& outz,
                   double& outw) const;

    /** This function applies the current transform to a position vector specified in regular
        coordinates and returns the transformed position vector in regular coordinates. If the
        transformed position is at infinity, the value of the returned vector is undefined. */
    Vec transform(Vec p) const;

private:
    double M[4][4];
};

//////////////////////////////////////////////////////////////////////

#endif

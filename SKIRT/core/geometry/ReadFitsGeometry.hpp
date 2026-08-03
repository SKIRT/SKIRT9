/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef READFITSGEOMETRY_HPP
#define READFITSGEOMETRY_HPP

#include "Array.hpp"
#include "GenGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The ReadFitsGeometry class is a subclass of the GenGeometry class, and describes a geometry
    characterized by observations. A 2D observed image FITS file is read into <tt>SKIRT</tt>, and
    deprojected assuming a certain position angle and inclination. The density is assumed
    to follow a exponential profile in the vertical directions, \f[ \rho(z) = \rho_0\,
    \exp\left(-\frac{|z|}{h_z}\right). \f] By running a <tt>SKIRT</tt> simulation with inclination
    of 0 degrees and position angle of the simulated galaxy, the <tt>SKIRT</tt> model images can be
    compared with the observations. The model geometry is set by nine parameters: the input
    filename, the pixel scale \f$pix\f$, the position angle \f$pa\f$, the inclination \f$incl\f$,
    the number of pixels in x and y direction \f$n_x\f$ and \f$n_y\f$, the center of galaxy in
    (x,y) image coordinates \f$x_c\f$ and \f$y_c\f$ and the vertical scale height \f$h_z\f$. */
class ReadFitsGeometry : public GenGeometry
{
    ITEM_CONCRETE(ReadFitsGeometry, GenGeometry, "a geometry read from a 2D FITS file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ReadFitsGeometry, "Level2")

        PROPERTY_STRING(filename, "the name of the input image file")

        PROPERTY_DOUBLE(pixelScale, "the physical scale of the image (length per pixel)")
        ATTRIBUTE_QUANTITY(pixelScale, "length")
        ATTRIBUTE_MIN_VALUE(pixelScale, "]0")

        PROPERTY_DOUBLE(positionAngle, "the position angle ω of the system")
        ATTRIBUTE_QUANTITY(positionAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(positionAngle, "-360 deg")
        ATTRIBUTE_MAX_VALUE(positionAngle, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(positionAngle, "0 deg")

        PROPERTY_DOUBLE(inclination, "the inclination angle θ of the system")
        ATTRIBUTE_QUANTITY(inclination, "posangle")
        ATTRIBUTE_MIN_VALUE(inclination, "[0 deg")
        ATTRIBUTE_MAX_VALUE(inclination, "90 deg[")
        ATTRIBUTE_DEFAULT_VALUE(inclination, "0 deg")

        PROPERTY_INT(numPixelsX, "number of pixels in the x direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "1")

        PROPERTY_DOUBLE(centerX, "x coordinate of the center (in pixels)")
        ATTRIBUTE_MIN_VALUE(centerX, "]0")

        PROPERTY_INT(numPixelsY, "number of pixels in the y direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "1")

        PROPERTY_DOUBLE(centerY, "y coordinate of the center (in pixels)")
        ATTRIBUTE_MIN_VALUE(centerY, "]0")

        PROPERTY_DOUBLE(scaleHeight, "the scale height")
        ATTRIBUTE_QUANTITY(scaleHeight, "length")
        ATTRIBUTE_MIN_VALUE(scaleHeight, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies the validity of the pixel scale, the inclination angle, the number of
        pixels in the x and y direction, the center of the image in x and y coordinates and the
        vertical scale height \f$h_z\f$. A vector of normalized cumulative pixel luminosities is
        computed, satisfying the condition that the total mass equals 1. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho(x,y,z)\f$ at the position (x,y,z). */
    double density(Position bfr) const override;

    /** This function generates a random position \f$(x,y,z)\f$ from the geometry, by drawing a random point
        from the appropriate probability density distribution function. The \f$(x,y)\f$ coordinates are
        derived from the normalized cumulative luminosity vector of the observed 2D projection. The z
        coordinate is derived from the vertical exponential probability distribution function. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density along the
        entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0) \, \text{d} x \f] */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density along the
        entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0) \, \text{d} y \f] */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z) \, \text{d} z \f] */
    double SigmaZ() const override;

private:
    /** This function transform a set of x and y coordinates in the inclined plane into the x and y
        coordinates in the (rotated) plane of the image. The transformation matrix \f$\bf{T}\f$ is
        given by:
        \f[
        {\bf{T}}
        =
        \begin{pmatrix}
        \cos\omega & -\sin\omega & 0 \\ \sin\omega & \cos\omega & 0 \\ 0 & 0 & 1
        \end{pmatrix}
        \begin{pmatrix}
        0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1
        \end{pmatrix}
        \f]
        where \f$\omega\f$ is the position angle.
    */
    void rotate(double& x, double& y) const;

    /** This function transforms a set of x and y coordinates in the plane of the image into the x and
        y coordinates in the (derotated) inclined plane. The transformation matrix is the inverse of
        the matrix used for the rotate function. */
    void derotate(double& x, double& y) const;

    /** This function projects the x coordinate of a point in the world coordinate system into the x
        coordinate of that point in the inclined plane. The following operation is applied:
        \f[x'=x\cdot \cos{\theta}\f]
        where \f$\theta\f$ is the inclination. */
    void project(double& x) const;

    /** This function deprojects the x coordinate of a point in the inclined plane into the x
        coordinate of that point in the world coordinate system. The following operation is applied:
        \f[x=\frac{x'}{\cos{\theta}}\f]
        where \f$\theta\f$ is the inclination. */
    void deproject(double& x) const;

    //======================== Data Members ========================

private:
    // the input image and the cumulative distribution
    Array _image;
    int _nx{0}, _ny{0}, _nz{0};
    Array _Xv;

    // other data members initialized during setup
    double _xmax{0.}, _ymax{0.}, _xmin{0.}, _ymin{0.};
    double _cospa{0.}, _sinpa{0.}, _cosi{0.}, _sini{0.};
    double _deltax{0.};
    double _C1x{0.}, _C1y{0.}, _C2x{0.}, _C2y{0.}, _C3x{0.}, _C3y{0.}, _C4x{0.}, _C4y{0.};
};

////////////////////////////////////////////////////////////////////

#endif

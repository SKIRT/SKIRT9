/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ROTATEGEOMETRYDECORATOR_HPP
#define ROTATEGEOMETRYDECORATOR_HPP

#include "GenGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** The RotateGeometryDecorator class is a decorator that applies an arbitrary rotation to any
    geometry. For the rotation, we use the general framework of the three Euler angles that can
    be used to decompose any rotation into a sequence of three individual rotations over the
    principle axes. We apply the following set of rotations (the so-called X-convention):
    <ul>
    <li> the first rotation is by an angle \f$\alpha\f$ about the Z axis.
    <li> the second rotation is by an angle \f$\beta\f$ about the new X' axis.
    <li> the third rotation is by an angle \f$\gamma\f$ about the new Z'' axis.
    </ul>
    If the original position of a vector is denoted as \f${\bf{r}}_{\text{orig}}\f$, the
    new position can be found as \f${\bf{r}} = {\bf{R}}\,{\bf{r}}_{\text{orig}}\f$, where
    the rotation matrix \f${\bf{R}}\f$ is given by
    \f[
    {\bf{R}}
    =
    \begin{pmatrix}
    \cos\gamma & \sin\gamma & 0 \\ -\sin\gamma & \cos\gamma & 0 \\ 0 & 0 & 1
    \end{pmatrix}
    \begin{pmatrix}
    1 & 0 & 0 \\ 0 & \cos\beta & \sin\beta \\ 0 & -\sin\beta & \cos\beta
    \end{pmatrix}
    \begin{pmatrix}
    \cos\alpha & \sin\alpha & 0 \\ -\sin\alpha & \cos\alpha & 0 \\ 0 & 0 & 1
    \end{pmatrix}
    \f]
    or explicitly
    \f[
    {\bf{R}}
    =
    \begin{pmatrix}
    \cos\alpha\cos\gamma-\sin\alpha\cos\beta\sin\gamma &
    \cos\gamma\sin\alpha+\cos\alpha\cos\beta\sin\gamma &
    \sin\beta\sin\gamma \\
    -\cos\alpha\sin\gamma-\sin\alpha\cos\beta\cos\gamma &
    -\sin\alpha\sin\gamma+\cos\alpha\cos\beta\cos\gamma &
    \sin\beta\cos\gamma \\
    \sin\alpha\sin\beta &
    -\cos\alpha\sin\beta &
    \cos\beta
    \end{pmatrix}
    \f]
    The properties of a RotateGeometryDecorator object are a reference to the Geometry object
    being decorated, and the three Euler angles \f$(\alpha,\beta,\gamma)\f$ that describe the
    rotation. The resulting geometry is identical to the geometry being decorated, except that
    the density distribution is rotated over the three Euler angles.
*/
class RotateGeometryDecorator : public GenGeometry
{
    ITEM_CONCRETE(RotateGeometryDecorator, GenGeometry, "a decorator that adds a rotation to any geometry")

        PROPERTY_ITEM(geometry, Geometry, "the geometry to be rotated")

        PROPERTY_DOUBLE(eulerAlpha, "the first Euler angle α")
        ATTRIBUTE_QUANTITY(eulerAlpha, "posangle")
        ATTRIBUTE_MIN_VALUE(eulerAlpha, "0 deg")
        ATTRIBUTE_MAX_VALUE(eulerAlpha, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(eulerAlpha, "0 deg")

        PROPERTY_DOUBLE(eulerBeta, "the second Euler angle β")
        ATTRIBUTE_QUANTITY(eulerBeta, "posangle")
        ATTRIBUTE_MIN_VALUE(eulerBeta, "0 deg")
        ATTRIBUTE_MAX_VALUE(eulerBeta, "180 deg")
        ATTRIBUTE_DEFAULT_VALUE(eulerBeta, "0 deg")

        PROPERTY_DOUBLE(eulerGamma, "the third Euler angle γ")
        ATTRIBUTE_QUANTITY(eulerGamma, "posangle")
        ATTRIBUTE_MIN_VALUE(eulerGamma, "0 deg")
        ATTRIBUTE_MAX_VALUE(eulerGamma, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(eulerGamma, "0 deg")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function calculates and stores some auxiliary values. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. It calls the density() function for the geometry being decorated with
        the derotated position \f${\bf{r}}_{\text{orig}} = {\bf{R}}^{\text{T}}\,{\bf{r}}\f$ as
        the argument. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. It calls the density() function for the geometry
        being decorated and rotates the resulting position \f${\bf{r}}_{\text{orig}}\f$. */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density
        along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f]
        It is impossible to calculate this value for a random value of the rotation angles. We simply
        return the X-axis surface density of the original geometry. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density
        along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f]
        It is impossible to calculate this value for a random value of the rotation angles. We simply
        return the Y-axis surface density of the original geometry. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density
        along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f]
        It is impossible to calculate this value for a random value of the rotation angles. We simply
        return the Z-axis surface density of the original geometry. */
    double SigmaZ() const override;

private:
    /** This function rotates a position \f${\bf{r}}_{\text{orig}}\f$, i.e.\ it returns the rotated
        position vector \f${\bf{r}} = {\bf{R}}\,{\bf{r}}_{\text{orig}}\f$. */
    Position rotate(Position bfrorig) const;

    /** This function derotates a position \f${\bf{r}}\f$, i.e.\ it returns the derotated
        position vector \f${\bf{r}}_{\text{orig}} = {\bf{R}}^{\text{T}}\,{\bf{r}}\f$. */
    Position derotate(Position bfr) const;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _sinalpha{0.}, _cosalpha{0.}, _sinbeta{0.}, _cosbeta{0.}, _singamma{0.}, _cosgamma{0.};
    double _R11{0.}, _R12{0.}, _R13{0.}, _R21{0.}, _R22{0.}, _R23{0.}, _R31{0.}, _R32{0.}, _R33{0.};
};

////////////////////////////////////////////////////////////////////

#endif

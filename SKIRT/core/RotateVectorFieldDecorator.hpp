/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ROTATEVECTORFIELDDECORATOR_HPP
#define ROTATEVECTORFIELDDECORATOR_HPP

#include "VectorField.hpp"

////////////////////////////////////////////////////////////////////

/** The RotateVectorFieldDecorator class is a decorator that applies an arbitrary rotation to any
    vector field. For the rotation, we use the general framework of the three Euler angles that can
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
    The properties of a RotateVectorFieldDecorator object are a reference to the vector field
    object being decorated, and the three Euler angles \f$(\alpha,\beta,\gamma)\f$ that describe
    the rotation. The resulting vector field is identical to the vector field being decorated,
    except that it is rotated over the three Euler angles. */
class RotateVectorFieldDecorator : public VectorField
{
    ITEM_CONCRETE(RotateVectorFieldDecorator, VectorField, "a decorator that adds a rotation to any vector field")
        ATTRIBUTE_TYPE_INSERT(RotateVectorFieldDecorator, "Dimension3")

        PROPERTY_ITEM(vectorField, VectorField, "the vector field to be rotated")

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
    /** This function returns the dimension of the vector field, which is 3 for this class because
        there are no guaranteed symmetries. */
    int dimension() const override;

    /** This function returns the value of the vector field at the position \f${\bf{r}}\f$. It
        calls the vector() function for the vector field being decorated with
        the derotated position \f${\bf{r}}_{\text{orig}} = {\bf{R}}^{\text{T}}\,{\bf{r}}\f$ as
        the argument. */
    Vec vector(Position bfr) const override;

private:
    /** This function rotates a vector \f${\bf{r}}_{\text{orig}}\f$, i.e.\ it returns the rotated
        vector \f${\bf{r}} = {\bf{R}}\,{\bf{r}}_{\text{orig}}\f$. */
    Vec rotate(Vec bfrorig) const;

    /** This function derotates a position \f${\bf{r}}\f$, i.e.\ it returns the derotated
        position \f${\bf{r}}_{\text{orig}} = {\bf{R}}^{\text{T}}\,{\bf{r}}\f$. */
    Position derotate(Position bfr) const;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _R11{0.}, _R12{0.}, _R13{0.}, _R21{0.}, _R22{0.}, _R23{0.}, _R31{0.}, _R32{0.}, _R33{0.};
};

////////////////////////////////////////////////////////////////////

#endif

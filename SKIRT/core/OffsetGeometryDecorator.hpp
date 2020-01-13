/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OFFSETGEOMETRYDECORATOR_HPP
#define OFFSETGEOMETRYDECORATOR_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The OffsetGeometryDecorator class is a decorator that adds an arbitrary offset to any geometry.
    The properties of an OffsetGeometryDecorator object include (1) a reference to the Geometry
    object being decorated and (2) three offsets in the x, y, and z directions. The resulting
    geometry is identical to the geometry being decorated, except that the density distribution is
    shifted over the specified offset. The geometry implemented by an OffsetGeometryDecorator
    object is 2D (axial symmetry) or 3D (no symmetries) depending on the symmetries of the geometry
    being decorated and on the specified offset. Specifically, it is 2D if the geometry being
    decorated is 1D or 2D and the offsets in the x and y directions are both zero. It is 3D if the
    geometry being decorated is 3D, or if at least one of the offsets in the x and y directions is
    nonzero. */
class OffsetGeometryDecorator : public Geometry
{
    ITEM_CONCRETE(OffsetGeometryDecorator, Geometry, "a decorator that adds an offset to any geometry")

        PROPERTY_ITEM(geometry, Geometry, "the geometry to be offset")

        PROPERTY_DOUBLE(offsetX, "the offset in the x direction")
        ATTRIBUTE_QUANTITY(offsetX, "length")
        ATTRIBUTE_DEFAULT_VALUE(offsetX, "0")
        ATTRIBUTE_INSERT(offsetX, "offsetX:Dimension3")

        PROPERTY_DOUBLE(offsetY, "the offset in the y direction")
        ATTRIBUTE_QUANTITY(offsetY, "length")
        ATTRIBUTE_DEFAULT_VALUE(offsetY, "0")
        ATTRIBUTE_INSERT(offsetY, "offsetY:Dimension3")

        PROPERTY_DOUBLE(offsetZ, "the offset in the z direction")
        ATTRIBUTE_QUANTITY(offsetZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(offsetZ, "0")
        ATTRIBUTE_INSERT(offsetZ, "offsetZ:Dimension2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry. It is 2 if the dimension of the
        geometry being decorated is 1 or 2 and the offsets in the x and y directions are both zero.
        It is 3 if the dimension of the geometry being decorated is 3, or if at least one of the
        offsets in the x and y directions is nonzero. */
    int dimension() const override;

    /** This function returns the density \f$\rho({\bf{r}})\f$ at the position
        \f${\bf{r}}\f$. It calls the density() function for the geometry being decorated with the
        translated position \f${\bf{r}}-{\bf{r}_\text{offset}}\f$. */
    double density(Position bfr) const override;

    /** This function generates a random position from the geometry, by drawing a random
        point from the three-dimensional probability density \f$p({\bf{r}})\, {\text{d}}{\bf{r}} =
        \rho({\bf{r}})\, {\text{d}}{\bf{r}}\f$. It calls the density() function for the geometry
        being decorated and returns the translated position \f${\bf{r}}+{\bf{r}_\text{offset}}\f$.
        */
    Position generatePosition() const override;

    /** This function returns the X-axis surface density, i.e. the integration of the density
        along the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f]
        It is impossible to calculate this value for a random value of the offset. The best
        option we have is to return the X-axis surface density of the original geometry, which
        is the true value in case there is only an offset in the x direction. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density
        along the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f]
        It is impossible to calculate this value for a random value of the offset. The best
        option we have is to return the Y-axis surface density of the original geometry, which
        is the true value in case there is only an offset in the y direction. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density
        along the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f]
        It is impossible to calculate this value for a random value of the offset. The best
        option we have is to return the Z-axis surface density of the original geometry, which
        is the true value in case there is only an offset in the z direction. */
    double SigmaZ() const override;
};

////////////////////////////////////////////////////////////////////

#endif

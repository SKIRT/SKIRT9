/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OFFSETVECTORFIELDDECORATOR_HPP
#define OFFSETVECTORFIELDDECORATOR_HPP

#include "VectorField.hpp"

////////////////////////////////////////////////////////////////////

/** The OffsetVectorFieldDecorator class is a decorator that adds an arbitrary offset to any vector
    field. The properties of an OffsetVectorFieldDecorator object include (1) a reference to the
    VectorField object being decorated and (2) three offsets in the x, y, and z directions. The
    resulting vector field is identical to the vector field being decorated, except that it is
    shifted over the specified offset. The vector field implemented by an
    OffsetVectorFieldDecorator object is 2D (axial symmetry) or 3D (no symmetries) depending on the
    symmetries of the vector field being decorated and on the specified offset. Specifically, it is
    2D if the vector field being decorated is 1D or 2D and the offsets in the x and y directions
    are both zero. It is 3D if the vector field being decorated is 3D, or if at least one of the
    offsets in the x and y directions is nonzero. */
class OffsetVectorFieldDecorator : public VectorField
{
    ITEM_CONCRETE(OffsetVectorFieldDecorator, VectorField, "a decorator that adds an offset to any vector field")

        PROPERTY_ITEM(vectorField, VectorField, "the vector field to be offset")

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
    /** This function returns the dimension of the vector field. It is 2 if the dimension of the
        vector field being decorated is 1 or 2 and the offsets in the x and y directions are both
        zero. It is 3 if the dimension of the vector field being decorated is 3, or if at least one
        of the offsets in the x and y directions is nonzero. */
    int dimension() const override;

    /** This function returns the value of the vector field at the position \f${\bf{r}}\f$. It
        calls the vector() function for the vector field being decorated with the translated
        position \f${\bf{r}}-{\bf{r}_\text{offset}}\f$. */
    Vec vector(Position bfr) const override;
};

////////////////////////////////////////////////////////////////////

#endif

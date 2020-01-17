/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDRICALVECTORFIELD_HPP
#define CYLINDRICALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** CylindricalVectorField represents a vector field with unit vectors rotating around the z-axis
    in a plane parallel to the xy-plane. The orientation is clockwise when looking along with the
    positive z-axis, and counterclockwise when looking down from the positive z-axis towards the
    xy-plane. */
class CylindricalVectorField : public VectorField
{
    ITEM_CONCRETE(CylindricalVectorField, VectorField, "a vector field rotating clockwise around the z-axis")
        ATTRIBUTE_TYPE_INSERT(CylindricalVectorField, "Dimension3")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the vector field, which is 3 for this class,
        indicating no symmetries (the vectors point in a different direction at each position). */
    int dimension() const override;

    /** This function returns a unit vector that is parallel with the xy-plane and orthogonal to
        the cylindrical radius component of the given position \f$\bf{r}\f$ with counterclockwise
        orientation in the xy-plane. Specifically, if \f$\bf{r}=(x,y,z)\f$, the function returns a
        normalized version of \f$(-y,x,0)\f$ unless \f$x=y=0\f$ in which case it returns the null
        vector. */
    Vec vector(Position bfr) const override;
};

////////////////////////////////////////////////////////////////////

#endif

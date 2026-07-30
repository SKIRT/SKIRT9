/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef UNIDIRECTIONALVECTORFIELD_HPP
#define UNIDIRECTIONALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** UnidirectionalVectorField represents a vector field with unit vectors that all point in the
    same direction. This direction is configured by the user; the components are automatically
    normalized to form a unit vector if needed. */
class UnidirectionalVectorField : public VectorField
{
    ITEM_CONCRETE(UnidirectionalVectorField, VectorField, "a vector field uniformly pointing in a given direction")

        PROPERTY_DOUBLE(fieldX, "the field direction, x component")
        ATTRIBUTE_DEFAULT_VALUE(fieldX, "0")
        ATTRIBUTE_INSERT(fieldX, "InMedia&fieldX:MediaDimension3,Dimension3;fieldX:SourceDimension3,Dimension3")

        PROPERTY_DOUBLE(fieldY, "the field direction, y component")
        ATTRIBUTE_DEFAULT_VALUE(fieldY, "0")
        ATTRIBUTE_INSERT(fieldY, "InMedia&fieldY:MediaDimension3,Dimension3;fieldY:SourceDimension3,Dimension3")

        PROPERTY_DOUBLE(fieldZ, "the field direction, z component")
        ATTRIBUTE_DEFAULT_VALUE(fieldZ, "1")
        ATTRIBUTE_INSERT(fieldZ, "InMedia&fieldZ:MediaDimension2,Dimension2;fieldZ:SourceDimension2,Dimension2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the configured direction is not the null vector, and it stores
        the normalized direction vector for later use. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the vector field, which for this class is 2 if the
        configured direction is parallel to the z-axis (indicating axial symmetry), and 3 otherwise
        (no symmetry). */
    int dimension() const override;

    /** This function returns a unit vector corresponding to the direction configured by the user.
        */
    Vec vector(Position bfr) const override;

    //======================== Data Members ========================

private:
    Vec _d;
};

////////////////////////////////////////////////////////////////////

#endif

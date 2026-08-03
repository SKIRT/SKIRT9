/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MESHGEOMETRY_HPP
#define MESHGEOMETRY_HPP

#include "Box.hpp"
#include "ImportedGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** MeshGeometry is an abstract class for representing a 3D geometry with a spatial density
    distribution that is discretized on some structured or unstructured tessellation of a cuboidal
    spatial domain. The class derives from the ImportedGeometry class, and its main (or only)
    function is to allow the user to configure the extent of the cuboidal domain of the
    tessellation, and to indicate the type of mass quantity to be imported for each cell.
    Subclasses need to define the actual tessellation being used, and deal with the other
    requirements set by the ImportedGeometry class. */
class MeshGeometry : public ImportedGeometry
{
    /** The enumeration type indicating the type of mass quantity to be imported. */
    ENUM_DEF(MassType, MassDensity, Mass, NumberDensity, Number)
        ENUM_VAL(MassType, MassDensity, "mass density")
        ENUM_VAL(MassType, Mass, "mass (volume-integrated mass density)")
        ENUM_VAL(MassType, NumberDensity, "number density")
        ENUM_VAL(MassType, Number, "number (volume-integrated number density)")
    ENUM_END()

    ITEM_ABSTRACT(MeshGeometry, ImportedGeometry, "a geometry imported from mesh-based data")

        PROPERTY_DOUBLE(minX, "the start point of the domain in the X direction")
        ATTRIBUTE_QUANTITY(minX, "length")

        PROPERTY_DOUBLE(maxX, "the end point of the domain in the X direction")
        ATTRIBUTE_QUANTITY(maxX, "length")

        PROPERTY_DOUBLE(minY, "the start point of the domain in the Y direction")
        ATTRIBUTE_QUANTITY(minY, "length")

        PROPERTY_DOUBLE(maxY, "the end point of the domain in the Y direction")
        ATTRIBUTE_QUANTITY(maxY, "length")

        PROPERTY_DOUBLE(minZ, "the start point of the domain in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

        PROPERTY_DOUBLE(maxZ, "the end point of the domain in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

        PROPERTY_ENUM(massType, MassType, "the type of mass quantity to be imported")
        ATTRIBUTE_DEFAULT_VALUE(massType, "MassDensity")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the configured domain has a positive volume. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

protected:
    /** This function returns the tessellation domain configured for this geometry. */
    const Box& domain() const { return _domain; }

    //======================== Data Members ========================

private:
    Box _domain;
};

////////////////////////////////////////////////////////////////////

#endif

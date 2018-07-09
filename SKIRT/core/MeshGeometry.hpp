/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MESHGEOMETRY_HPP
#define MESHGEOMETRY_HPP

#include "ImportedGeometry.hpp"
#include "Box.hpp"

////////////////////////////////////////////////////////////////////

/** MeshGeometry is an abstract class for representing a 3D geometry with a spatial density
    distribution that is discretized on some structured or unstructured tessellation of a cuboidal
    spatial domain. The class derives from the ImportedGeometry class, and its main (or only)
    function is to allow the user to configure the extent of the cuboidal domain of the
    tessellation, and to indicate whether the mass or the density is being specified for each cell.
    Subclasses need to define the actual tessellation being used, and deal with the other
    requirements set by the ImportedGeometry class. */
class MeshGeometry : public ImportedGeometry
{
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

    PROPERTY_BOOL(useMass, "import the integrated mass for each cell (rather than the mass density)")
        ATTRIBUTE_DEFAULT_VALUE(useMass, "false")

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

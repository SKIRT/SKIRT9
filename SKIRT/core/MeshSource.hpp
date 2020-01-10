/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MESHSOURCE_HPP
#define MESHSOURCE_HPP

#include "Box.hpp"
#include "ImportedSource.hpp"

////////////////////////////////////////////////////////////////////

/** MeshSource is an abstract class for representing a primary radiation source with a luminosity
    distribution that is discretized on some structured or unstructured tessellation of a cuboidal
    spatial domain. The class derives from the ImportedSource class, and its main (or only)
    function is to allow the user to configure the extent of the cuboidal domain of the
    tessellation. Subclasses need to define the actual tessellation being used, and deal with the
    other requirements set by the ImportedSource class. */
class MeshSource : public ImportedSource
{
    ITEM_ABSTRACT(MeshSource, ImportedSource, "a primary source imported from mesh-based data")

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

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the configured domain has a positive volume. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

protected:
    /** This function returns the tessellation domain configured for this source. */
    const Box& domain() const { return _domain; }

    //======================== Data Members ========================

private:
    Box _domain;
};

////////////////////////////////////////////////////////////////////

#endif

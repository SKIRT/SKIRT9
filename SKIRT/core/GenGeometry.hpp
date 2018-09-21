/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GENGEOMETRY_HPP
#define GENGEOMETRY_HPP

#include "Geometry.hpp"

////////////////////////////////////////////////////////////////////

/** The GenGeometry class is an abstract subclass of the Geometry class, and serves
    as a base class for general geometries without a specific symmetry. */
class GenGeometry : public Geometry
{
    ITEM_ABSTRACT(GenGeometry, Geometry, "a geometry without a specific symmetry")
        ATTRIBUTE_TYPE_INSERT(GenGeometry, "Dimension3")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry, which is 3 for all subclasses of this
        class since it is a base class for geometries without a specific symmetry. */
    int dimension() const override;
};

////////////////////////////////////////////////////////////////////

#endif

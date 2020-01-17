/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VECTORFIELD_HPP
#define VECTORFIELD_HPP

#include "Position.hpp"
#include "SimulationItem.hpp"

//////////////////////////////////////////////////////////////////////

/** VectorField is an abstract base class for describing vector fields. Specifically, an instance
    of a VectorField subclass describes a normalized vector field, i.e. a spatial distribution of
    3D vector values. The field is normalized so that the maximum length (norm) of the vectors in
    the field is equal to one. This also implies that the vector components are dimensionless. As a
    result, vector fields can be used for various purposes, including the specification of magnetic
    fields and velocity fields.

    This base class defines the interface for all subclasses in the hierarchy. It consists of the
    key vector(Position) function, which returns the value of the vector field at a given position,
    and the dimension() function, which returns a value indicating the (lack of) spatial symmetry
    in the field. Each concrete subclass must implement these functions. */
class VectorField : public SimulationItem
{
    ITEM_ABSTRACT(VectorField, SimulationItem, "a vector field")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the vector field, which depends on its (lack of)
        symmetry. A value of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of
        these symmetries. The function's implementation must be provided in a subclass. */
    virtual int dimension() const = 0;

    /** This function returns the value of the vector field (a 3D vector with dimensionless
        components) at the position \f${\bf{r}}\f$. The function's implementation must be provided
        in a subclass. */
    virtual Vec vector(Position bfr) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

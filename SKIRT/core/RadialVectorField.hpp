/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIALVECTORFIELD_HPP
#define RADIALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** RadialVectorField represents a vector field with unit vectors pointing away from the origin. */
class RadialVectorField : public VectorField
{
    ITEM_CONCRETE(RadialVectorField, VectorField, "a vector field pointing away from the origin")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the vector field, which is 1 for this class,
        indicating spherical symmetry. */
    int dimension() const override;

    /** This function returns a unit vector pointing away from the origin at the given position
        \f$\bf{r}\f$. Specifically, it returns \f$\bf{r}/||\bf{r}||\f$ unless \f$||\bf{r}||=0\f$ in
        which case it returns the null vector. */
    Vec vector(Position bfr) const override;
};

////////////////////////////////////////////////////////////////////

#endif

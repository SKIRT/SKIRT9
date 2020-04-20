/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIALVECTORFIELD_HPP
#define RADIALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** RadialVectorField represents a vector field with vectors pointing away from the origin with
    lengths that depend on the spatial location and on a scale radius \f$r_1\f$ specified by the
    user as a property. Specifically, the magnitude \f$v\f$ of the vectors is given by \f[ v =
    \begin{cases} r/r_1 & \mathrm{for}\;0\le r<r_1 \\ 1 & \mathrm{for}\;r\ge r_1 \end{cases} \f] or
    by \f$v=1\f$ if \f$r_1=0\f$ (the default value). In other words, for the default transition
    radius of zero, the field consists of unit vectors. For other transition radii, the vector
    magnitude scales from zero at \f$r=0\f$ to unity at \f$r=r_1\f$, and stays at unity for
    \f$r>r_1\f$. */
class RadialVectorField : public VectorField
{
    ITEM_CONCRETE(RadialVectorField, VectorField, "a vector field pointing away from the origin")
        ATTRIBUTE_TYPE_INSERT(RadialVectorField, "Dimension3")

        PROPERTY_DOUBLE(unityRadius, "the radius where the magnitude of the vectors is unity")
        ATTRIBUTE_QUANTITY(unityRadius, "length")
        ATTRIBUTE_MIN_VALUE(unityRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(unityRadius, "0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the vector field, which is 3 for this class,
        indicating no symmetries (the vectors point in a different direction at each position). */
    int dimension() const override;

    /** This function returns a unit vector pointing away from the origin at the given position
        \f$\bf{r}\f$. Specifically, it returns \f$\bf{r}/||\bf{r}||\f$ unless \f$||\bf{r}||=0\f$ in
        which case it returns the null vector. */
    Vec vector(Position bfr) const override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIALVECTORFIELD_HPP
#define RADIALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** RadialVectorField represents a vector field with vectors pointing away from the origin and with
    a magnitude that varies with radial distance from the origin. The class has two configuration
    properties: the distance at which the vector magnitude becomes unity and a power law exponent.
    If the exponent is positive, the magnitude ramps up from zero at the origin to unity at the
    specified distance, and remains unity from there on. If the power law index is negative, the
    magnitude is unity from the origin to the specified distance, and then ramps down from there to
    zero at infinity.

    Given the user-configurable unity radius \f$r_1\ge 0\f$ and power-law exponent \f$\alpha\f$,
    the magnitude \f$v(r)\f$ of the vectors as a function of radial distance \f$r\f$ is given by

    \f[ r_1=0 \lor \alpha=0 \qquad v(r) = \begin{cases} 0 & \mathrm{for}\;r=0 \\ 1 &
    \mathrm{for}\;r>0 \end{cases} \f]

    \f[ r_1>0 \land \alpha>0 \qquad v(r) = \begin{cases} 0 & \mathrm{for}\;r=0 \\ (r/r_1)^\alpha &
    \mathrm{for}\;0<r<r_1 \\ 1 & \mathrm{for}\;r\ge r_1 \end{cases} \f]

    \f[ r_1>0 \land \alpha<0 \qquad v(r) = \begin{cases} 0 & \mathrm{for}\;r=0 \\ 1 &
    \mathrm{for}\;0<r<r_1 \\ (r/r_1)^\alpha & \mathrm{for}\;r \ge r_1 \end{cases} \f]

    The first equation defines the vector field for the special values \f$r_1=0\f$ and/or
    \f$\alpha=0\f$ to consist solely of unit vectors except at the origin, where the direction of
    the vector is undefined. The default value of \f$r_1=0\f$ causes this behavior. If the unity
    radius is set to a value \f$r_1>0\f$, the default value of \f$\alpha=1\f$ specifies linear
    behavior: the vector magnitude scales linearly from zero at the origin to unity at \f$r=r_1\f$,
    and stays at unity for larger radii. */
class RadialVectorField : public VectorField
{
    ITEM_CONCRETE(RadialVectorField, VectorField, "a vector field pointing away from the origin")
        ATTRIBUTE_TYPE_INSERT(RadialVectorField, "Dimension3")

        PROPERTY_DOUBLE(unityRadius, "the radius where the magnitude of the vectors is unity")
        ATTRIBUTE_QUANTITY(unityRadius, "length")
        ATTRIBUTE_MIN_VALUE(unityRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(unityRadius, "0")

        PROPERTY_DOUBLE(exponent, "the power-law exponent governing the radial dependence of the magnitude")
        ATTRIBUTE_DEFAULT_VALUE(exponent, "1")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the vector field, which is 3 for this class,
        indicating no symmetries (the vectors point in a different direction at each position). */
    int dimension() const override;

    /** This function returns a vector pointing away from the origin at the given position
        \f$\bf{r}\f$ according to the definitions given in the class header. */
    Vec vector(Position bfr) const override;
};

////////////////////////////////////////////////////////////////////

#endif

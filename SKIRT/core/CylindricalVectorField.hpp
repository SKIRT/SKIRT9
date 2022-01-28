/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDRICALVECTORFIELD_HPP
#define CYLINDRICALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** CylindricalVectorField represents a vector field with vectors rotating around the z-axis in a
    plane parallel to the xy-plane. The orientation is clockwise when looking along with the
    positive z-axis, and counterclockwise when looking down from the positive z-axis towards the
    xy-plane.

    The magnitude of the vectors varies with radial distance from the z-axis. The class has two
    configuration properties: the distance at which the vector magnitude becomes unity and a power
    law exponent. If the exponent is positive, the magnitude ramps up from zero at the z-axis to
    unity at the specified distance, and remains unity from there on. If the power law index is
    negative, the magnitude is unity from the z-axis to the specified distance, and then ramps down
    from there to zero at infinity.

    Given the user-configurable unity radius \f$R_1\ge 0\f$ and power-law exponent \f$\alpha\f$,
    the magnitude \f$v(R)\f$ of the vectors as a function of radial distance \f$R=\sqrt{x^2+y^2}\f$
    is given by

    \f[ R_1=0 \lor \alpha=0 \qquad v(R) = \begin{cases} 0 & \mathrm{for}\;R=0 \\ 1 &
    \mathrm{for}\;R>0 \end{cases} \f]

    \f[ R_1>0 \land \alpha>0 \qquad v(R) = \begin{cases} 0 & \mathrm{for}\;R=0 \\ (R/R_1)^\alpha &
    \mathrm{for}\;0<R<R_1 \\ 1 & \mathrm{for}\;R\ge R_1 \end{cases} \f]

    \f[ R_1>0 \land \alpha<0 \qquad v(R) = \begin{cases} 0 & \mathrm{for}\;R=0 \\ 1 &
    \mathrm{for}\;0<R<R_1 \\ (R/R_1)^\alpha & \mathrm{for}\;R \ge R_1 \end{cases} \f]

    The first equation defines the vector field for the special values \f$R_1=0\f$ and/or
    \f$\alpha=0\f$ to consist solely of unit vectors except at the z-axis, where the direction of
    the vector is undefined. The default value of \f$R_1=0\f$ causes this behavior. If the unity
    radius is set to a value \f$R_1>0\f$, the default value of \f$\alpha=1\f$ specifies linear
    behavior: the vector magnitude scales linearly from zero at the z-axis to unity at \f$R=R_1\f$,
    and stays at unity for larger radii. */
class CylindricalVectorField : public VectorField
{
    ITEM_CONCRETE(CylindricalVectorField, VectorField, "a vector field rotating clockwise around the z-axis")
        ATTRIBUTE_TYPE_INSERT(CylindricalVectorField, "Dimension3")

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

    /** This function returns a vector that is parallel with the xy-plane and orthogonal to the
        cylindrical radius component of the given position \f$\bf{r}\f$ according to the
        definitions given in the class header. */
    Vec vector(Position bfr) const override;
};

////////////////////////////////////////////////////////////////////

#endif

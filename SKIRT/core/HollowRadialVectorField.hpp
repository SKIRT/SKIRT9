/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HOLLOWRADIALVECTORFIELD_HPP
#define HOLLOWRADIALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** HollowRadialVectorField represents a vector field with vectors pointing away from the origin
    and with a magnitude that varies with radial distance from the origin. The vector magnitudes
    increase from zero at a given radius to unity at infinity. Within the given radius, all vector
    magnitudes are zero, creating a magnitude cavity to which this class owes its name.

    The class has two configuration properties: the distance \f$r_0\ge 0\f$ at which the vector
    magnitude starts increasing and a power law exponent \f$\alpha\ge 0\f$. Both properties must be
    non-negative. Given these properties, the magnitude \f$v(r)\f$ of the vectors as a function of
    radial distance \f$r\f$ is given by

    \f[ v(r) = \begin{cases} 0 & \mathrm{for}\;r\le r_0 \\ \left(1-\frac{r_0}{r}\right)^\alpha &
    \mathrm{for}\;r>r_0 \end{cases} \f]

    With the default value of \f$r_0=0\f$, there is no central cavity. The default value of
    \f$\alpha=0\f$ specifies a step function from zero to unity at \f$r=r_0\f$.

    This radial magnitude dependence is inspired by Braibant et al. 2017 (A&A, 607, A32) and Murray
    et al. 1995 (ApJ, 451, 498). It is intended to serve as a velocity field for model geometries
    that have no mass within a given radius, such as the dusty torus of an active galactic nucleus.
    */
class HollowRadialVectorField : public VectorField
{
    ITEM_CONCRETE(HollowRadialVectorField, VectorField,
                  "a vector field pointing away from the origin with a central cavity")
        ATTRIBUTE_TYPE_INSERT(HollowRadialVectorField, "Dimension3")

        PROPERTY_DOUBLE(zeroRadius, "the radius within which the magnitude of the vectors is zero")
        ATTRIBUTE_QUANTITY(zeroRadius, "length")
        ATTRIBUTE_MIN_VALUE(zeroRadius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(zeroRadius, "0")

        PROPERTY_DOUBLE(exponent, "the power-law exponent governing the radial dependence of the magnitude")
        ATTRIBUTE_MIN_VALUE(exponent, "[0")
        ATTRIBUTE_DEFAULT_VALUE(exponent, "0")

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

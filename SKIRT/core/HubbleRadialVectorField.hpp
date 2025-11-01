/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HUBBLERADIALVECTORFIELD_HPP
#define HUBBLERADIALVECTORFIELD_HPP

#include "VectorField.hpp"

//////////////////////////////////////////////////////////////////////

/** HubbleRadialVectorField represents a vector field with vectors pointing away from the origin
    and a magnitude that varies with radial distance from the origin. Up to the turnover radius
    \f$(r_t \ge 0)\f$, the velocity is constantly accelerating. Beyond the turnover radius, the
    velocity is constantly decelarating. Beyond the maximum radius \f$(r_\mathrm{max} \ge r_t)\f$,
    the velocity is zero:
 
    \f[ v(r) = \begin{cases}\frac{r}{r_t} & \mathrm{for}\;r\le r_t \\ \left(1-\frac{r - r_t}{r_\mathrm{max}-r_t}\right)
    & \mathrm{for}\;r_t < r\le r_\mathrm{max} \\0 & \mathrm{for}\;r> r_\mathrm{max} \end{cases} \f]
 
    This radial magnitude dependence was presented in Das et al. 2006 (AJ, 132, 620D) as the
    velocity field for the biconical polar outflow in NGC1068.
    */
class HubbleRadialVectorField : public VectorField
{
    ITEM_CONCRETE(HubbleRadialVectorField, VectorField,
                  "a Hubble flow pointing away from the origin with a constant acceleration and decelaration")
        ATTRIBUTE_TYPE_INSERT(HubbleRadialVectorField, "Dimension3")

        PROPERTY_DOUBLE(turnoverRadius, "the turnover radius marking the transition from acceleration to decelaration")
        ATTRIBUTE_QUANTITY(turnoverRadius, "length")
        ATTRIBUTE_MIN_VALUE(turnoverRadius, "[0")

        PROPERTY_DOUBLE(maxRadius, "the maximum radius beyond which the outflow velocity is zero")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "[0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that \f$r_t \le r_\mathrm{max}\f$. */
    void setupSelfBefore() override;

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

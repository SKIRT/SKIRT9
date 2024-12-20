
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHEPOWERLAWREDISTRIBUTEGEOMETRYDECORATOR_HPP
#define SPHEPOWERLAWREDISTRIBUTEGEOMETRYDECORATOR_HPP

#include "RedistributeGeometryDecorator.hpp"

////////////////////////////////////////////////////////////////////

/** The SphePowerLawRedistributeGeometryDecorator class implements a decorator that adjusts another
    geometry by multiplying the density with a spherical power law weight function \f[ \rho'(r,
    \theta, \phi) = n \rho(r, \theta, \phi) r^{-p}. \f] There is also a spherical clipping region
    around the origin determined by a radius \f$r_0\f$ where the density is made zero to cut out
    the singularity. */
class SphePowerLawRedistributeGeometryDecorator : public RedistributeGeometryDecorator
{
    ITEM_CONCRETE(SphePowerLawRedistributeGeometryDecorator, RedistributeGeometryDecorator,
                  "a decorator that redistributes another geometry with a spherical power law")

        PROPERTY_DOUBLE(exponent, "the negative power of the weight function")
        ATTRIBUTE_MIN_VALUE(exponent, "]0")

        PROPERTY_DOUBLE(minRadius, "the radius of the clipping sphere")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "]0 pc")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** The dimension of the geometry after applying the decorator cannot change and is thus the
        dimension of the original geometry. */
    int dimension() const override;

protected:
    /** The weight function is the power law: \f$r^{-p}\f$. */
    double weight(Position bfr) const override;

    /** The max weight, used in the rejection method, is equal to \f$r_0^{-p}\f$ with \f$r_0>0\f$
        the radius of the clipping sphere. */
    double maxWeight() const override;
};

////////////////////////////////////////////////////////////////////

#endif


/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef AXPOWERLAWREDISTRIBUTEGEOMETRYDECORATOR_HPP
#define AXPOWERLAWREDISTRIBUTEGEOMETRYDECORATOR_HPP

#include "RedistributeGeometryDecorator.hpp"

////////////////////////////////////////////////////////////////////

/** The AxPowerLawRedistributeGeometryDecorator class implements a decorator that adjusts another
    geometry by multiplying the density with an axial power law weight function \f[ \rho'(R, \phi,
    z) = n \rho(R, \phi, z) R^{-p_R}z^{-p_z}. \f] There are also two clipping regions, one around
    the z-axis determined by a radius \f$R_0\f$ and one above and below the xy-plane determined by
    \f$z=\pm z_0\f$, where the density is made zero to cut out the singularity. If an exponent is
    chosen to be zero, then the clipping region for that exponent will not be prompted during the
    Q&A. */
class AxPowerLawRedistributeGeometryDecorator : public RedistributeGeometryDecorator
{
    ITEM_CONCRETE(AxPowerLawRedistributeGeometryDecorator, RedistributeGeometryDecorator,
                  "a decorator that redistributes another geometry with an axial power law")
        ATTRIBUTE_TYPE_INSERT(AxPowerLawRedistributeGeometryDecorator, "Dimension2")

        PROPERTY_DOUBLE(RExponent, "the negative power of the radial part of the weight function")
        ATTRIBUTE_MIN_VALUE(RExponent, "[0")

        PROPERTY_DOUBLE(zExponent, "the negative power of the vertical part of the weight function")
        ATTRIBUTE_MIN_VALUE(zExponent, "[0")

        PROPERTY_DOUBLE(minR, "the radius of the clipping cylinder around the z-axis")
        ATTRIBUTE_QUANTITY(minR, "length")
        ATTRIBUTE_MIN_VALUE(minR, "[0 pc")
        ATTRIBUTE_RELEVANT_IF(minR, "RExponent")

        PROPERTY_DOUBLE(minZ, "half the height of the clipping region around the xy-plane")
        ATTRIBUTE_QUANTITY(minZ, "length")
        ATTRIBUTE_MIN_VALUE(minZ, "[0 pc")
        ATTRIBUTE_RELEVANT_IF(minZ, "zExponent")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function checks if not both the \f$R\f$ and \f$z\f$ exponents are 0. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** The dimension of the geometry after applying the decorator can only change spherical
        geometries into axially symmetric, otherwise it doesn't change the dimension. */
    int dimension() const override;

protected:
    /** The weight function is the power law: \f$R^{-p_R}z^{-p_z}\f$. */
    double weight(Position bfr) const override;

    /** The max weight is used in the rejection method and is equal to \f$R_0^{-p_R}z_0^{-p_z}\f$.
        */
    double maxWeight() const override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation and backwards compatibility
    const double& _R0{_minR};
    const double& _z0{_minZ};
    const double& _pR{_RExponent};
    const double& _pz{_zExponent};
};

////////////////////////////////////////////////////////////////////

#endif

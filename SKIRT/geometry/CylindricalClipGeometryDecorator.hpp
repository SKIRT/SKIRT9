/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDRICALCLIPGEOMETRYDECORATOR_HPP
#define CYLINDRICALCLIPGEOMETRYDECORATOR_HPP

#include "ClipGeometryDecorator.hpp"

////////////////////////////////////////////////////////////////////

/** The CylindricalClipGeometryDecorator class is a decorator that adjusts another geometry by
    setting the density equal to zero inside or outside an infinitely long cylinder centered at the
    origin and oriented along the Z-axis. The radius of the cylinder can be chosen. */
class CylindricalClipGeometryDecorator : public ClipGeometryDecorator
{
    ITEM_CONCRETE(CylindricalClipGeometryDecorator, ClipGeometryDecorator,
                  "a decorator that clips another geometry using a cylinder")
        ATTRIBUTE_TYPE_INSERT(CylindricalClipGeometryDecorator,
                              "InMedia:MediaDimension2,Dimension2;SourceDimension2,Dimension2")

        PROPERTY_DOUBLE(clipRadius, "the radius of the clipping cylinder")
        ATTRIBUTE_QUANTITY(clipRadius, "length")
        ATTRIBUTE_MIN_VALUE(clipRadius, "[0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry. If the original geometry has a dimension 3,
        so will the new geometry. Otherwise, i.e. if the original geometry is spherically or axisymmetric,
        the dimension is 2. */
    int dimension() const override;

    /** This function returns true if the specified position is inside the cylinder defined by the
        properties of this class. */
    bool inside(Position bfr) const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f]
        Because the clip region is defined purely in terms of cylindrical radius, the Z-axis itself
        (at R=0) is always either entirely inside or entirely outside the clipped region, so unlike
        the generic ClipGeometryDecorator::SigmaX() and ClipGeometryDecorator::SigmaY() this value
        can be computed exactly rather than approximated. If the inside region is being removed,
        the entire Z-axis lies in the removed region and this function returns zero; otherwise the
        entire Z-axis is retained and this function returns the corresponding value of the geometry
        being decorated, renormalized by the same factor used in density(). This class does not
        override SigmaX() or SigmaY(): a cylindrical clip partially clips the X-axis and Y-axis (at
        \f$|x|>{\text{clipRadius}}\f$ and \f$|y|>{\text{clipRadius}}\f$ respectively when removing
        the outside region), so no exact shortcut is available for those two, and the inherited
        ClipGeometryDecorator implementation is used instead. */
    double SigmaZ() const override;
};

////////////////////////////////////////////////////////////////////

#endif

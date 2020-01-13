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
        ATTRIBUTE_TYPE_INSERT(CylindricalClipGeometryDecorator, "Dimension2")

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

    /** This function returns the X-axis surface density, i.e. the integration of the density along
        the entire X-axis, \f[ \Sigma_X = \int_{-\infty}^\infty \rho(x,0,0)\,{\text{d}}x. \f] It
        returns the corresponding value of the geometry being decorated after normalization. */
    double SigmaX() const override;

    /** This function returns the Y-axis surface density, i.e. the integration of the density along
        the entire Y-axis, \f[ \Sigma_Y = \int_{-\infty}^\infty \rho(0,y,0)\,{\text{d}}y. \f] It
        returns the corresponding value of the geometry being decorated after normalization. */
    double SigmaY() const override;

    /** This function returns the Z-axis surface density, i.e. the integration of the density along
        the entire Z-axis, \f[ \Sigma_Z = \int_{-\infty}^\infty \rho(0,0,z)\,{\text{d}}z. \f] If
        the inside region is being removed, this function returns zero; otherwise it returns the
        corresponding value of the geometry being decorated. */
    double SigmaZ() const override;
};

////////////////////////////////////////////////////////////////////

#endif

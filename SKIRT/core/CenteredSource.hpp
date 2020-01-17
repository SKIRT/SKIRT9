/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CENTEREDSOURCE_HPP
#define CENTEREDSOURCE_HPP

#include "SpecialtySource.hpp"

//////////////////////////////////////////////////////////////////////

/** CenteredSource is an abstract class for representing specialty radiation sources that have a
    configurable central point in space, but cannot be implemented as a GeometricSource because of
    emission anisotropies. A typical example is a primary source that provides a homogeneous
    background radiation field in a well-defined spatial region.

    This class allows the user to configure the central location, and through its base class, it
    offer user configuration of the %SED and luminosity normalization, and possibly a bulk
    velocity. */
class CenteredSource : public SpecialtySource
{
    ITEM_ABSTRACT(CenteredSource, SpecialtySource, "a centered primary source")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CenteredSource, "Level2")

        PROPERTY_DOUBLE(centerX, "the center of the source, x component")
        ATTRIBUTE_QUANTITY(centerX, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerX, "0")
        ATTRIBUTE_INSERT(centerX, "centerX:Dimension3")

        PROPERTY_DOUBLE(centerY, "the center of the source, y component")
        ATTRIBUTE_QUANTITY(centerY, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerY, "0")
        ATTRIBUTE_INSERT(centerY, "centerY:Dimension3")

        PROPERTY_DOUBLE(centerZ, "the center of the source, z component")
        ATTRIBUTE_QUANTITY(centerZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(centerZ, "0")
        ATTRIBUTE_INSERT(centerZ, "centerZ:Dimension2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the geometry of the source, which is the same as the
        instrinsic dimension (to be provided by the subclass) except if the center is not in the
        coordinate system origin. */
    int geometryDimension() const override;

    //============== Functions to be implemented in each subclass =============

    /** This function returns the intrinsic dimension of the spatial distribution implemented by
        the subclass, taking into account anisotropic emission or polarization, if any, but
        assuming that the center is in the coordinate system origin (and ignoring any bulk
        velocity). The function must be implemented in a subclass. */
    virtual int intrinsicDimension() const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif

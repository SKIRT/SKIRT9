/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef AXISMATERIALNORMALIZATION_HPP
#define AXISMATERIALNORMALIZATION_HPP

#include "MaterialNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** AxisMaterialNormalization is an intermediate abstract class representing material normalization
    types that require the user to select one of the axes of the model coordinatye system. */
class AxisMaterialNormalization : public MaterialNormalization
{
    /** The enumeration type indicating the coordinate axis for normalization. */
    ENUM_DEF(Axis, X, Y, Z)
        ENUM_VAL(Axis, X, "the X axis of the model coordinate sytem")
        ENUM_VAL(Axis, Y, "the Y axis of the model coordinate sytem")
        ENUM_VAL(Axis, Z, "the Z axis of the model coordinate sytem")
    ENUM_END()

    ITEM_ABSTRACT(AxisMaterialNormalization, MaterialNormalization,
                  "normalization by defining some quantity along a coordinate axis")

        PROPERTY_ENUM(axis, Axis, "the axis along which to specify the normalization")
        ATTRIBUTE_DEFAULT_VALUE(axis, "Z")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** As service to subclasses, this function returns the column density of the specified
        geometry along the coordinate axis configured by the user. */
    double geometryColumnDensity(const Geometry* geom) const;
};

////////////////////////////////////////////////////////////////////

#endif

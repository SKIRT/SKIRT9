/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPTICALDEPTHMATERIALNORMALIZATION_HPP
#define OPTICALDEPTHMATERIALNORMALIZATION_HPP

#include "MaterialNormalization.hpp"

//////////////////////////////////////////////////////////////////////

/** A OpticalDepthMaterialNormalization object normalizes the amount of material in a geometric
    medium by specifying the optical depth along one of the coordinate axes at a given wavelength. */
class OpticalDepthMaterialNormalization : public MaterialNormalization
{
    /** The enumeration type indicating the coordinate axis for normalization. */
    ENUM_DEF(Axis, X, Y, Z)
    ENUM_VAL(Axis, X, "the X axis of the model coordinate sytem")
    ENUM_VAL(Axis, Y, "the Y axis of the model coordinate sytem")
    ENUM_VAL(Axis, Z, "the Z axis of the model coordinate sytem")
    ENUM_END()

    ITEM_CONCRETE(OpticalDepthMaterialNormalization, MaterialNormalization,
                  "normalization by defining the optical depth along a coordinate axis")

    PROPERTY_ENUM(axis, Axis, "the axis along which to specify the optical depth")
        ATTRIBUTE_DEFAULT_VALUE(axis, "Z")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    PROPERTY_DOUBLE(opticalDepth, "the optical depth along this axis at this wavelength")
    ATTRIBUTE_MIN_VALUE(opticalDepth, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the total number of entities and total mass in the medium, in that
        order, given a geometry and material mix in addition to the user configuration options offered
        by this class. */
    std::pair<double,double> numberAndMass(const Geometry* geom, const MaterialMix* mix) const override;
};

//////////////////////////////////////////////////////////////////////

#endif

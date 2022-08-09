/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPTICALDEPTHMATERIALNORMALIZATION_HPP
#define OPTICALDEPTHMATERIALNORMALIZATION_HPP

#include "AxisMaterialNormalization.hpp"
#include "MaterialWavelengthRangeInterface.hpp"

////////////////////////////////////////////////////////////////////

/** A OpticalDepthMaterialNormalization object normalizes the amount of material in a geometric
    medium by specifying the optical depth along one of the coordinate axes at a given wavelength.
    It implements the MaterialWavelengthRangeInterface to indicate that wavelength-dependent
    material properties will be required for the configured wavelength. */
class OpticalDepthMaterialNormalization : public AxisMaterialNormalization, public MaterialWavelengthRangeInterface
{
    ITEM_CONCRETE(OpticalDepthMaterialNormalization, AxisMaterialNormalization,
                  "normalization by defining the optical depth along a coordinate axis")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to specify the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

        PROPERTY_DOUBLE(opticalDepth, "the optical depth along this axis at this wavelength")
        ATTRIBUTE_MIN_VALUE(opticalDepth, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the total number of entities and total mass in the medium, in that
        order, given a geometry and material mix in addition to the user configuration options offered
        by this class. */
    std::pair<double, double> numberAndMass(const Geometry* geom, const MaterialMix* mix) const override;

    /** This function returns a wavelength range corresponding to the user-configured wavelength,
        indicating that wavelength-dependent material properties will be required for this
        wavelength. */
    Range wavelengthRange() const override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NUMBERCOLUMNMATERIALNORMALIZATION_HPP
#define NUMBERCOLUMNMATERIALNORMALIZATION_HPP

#include "AxisMaterialNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A NumberColumnMaterialNormalization object normalizes the amount of material in a geometric
    medium by specifying the number column density along one of the coordinate axes. */
class NumberColumnMaterialNormalization : public AxisMaterialNormalization
{
    ITEM_CONCRETE(NumberColumnMaterialNormalization, AxisMaterialNormalization,
                  "normalization by defining the number column density along a coordinate axis")

        PROPERTY_DOUBLE(numberColumnDensity, "the number column density along this axis")
        ATTRIBUTE_QUANTITY(numberColumnDensity, "numbersurfacedensity")
        ATTRIBUTE_MIN_VALUE(numberColumnDensity, "]0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the total number of entities and total mass in the medium, in that
        order, given a geometry and material mix in addition to the user configuration options offered
        by this class. */
    std::pair<double, double> numberAndMass(const Geometry* geom, const MaterialMix* mix) const override;
};

////////////////////////////////////////////////////////////////////

#endif

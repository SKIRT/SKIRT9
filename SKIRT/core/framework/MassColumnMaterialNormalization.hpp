/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MASSCOLUMNMATERIALNORMALIZATION_HPP
#define MASSCOLUMNMATERIALNORMALIZATION_HPP

#include "AxisMaterialNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A MassColumnMaterialNormalization object normalizes the amount of material in a geometric
    medium by specifying the mass column density along one of the coordinate axes. */
class MassColumnMaterialNormalization : public AxisMaterialNormalization
{
    ITEM_CONCRETE(MassColumnMaterialNormalization, AxisMaterialNormalization,
                  "normalization by defining the mass column density along a coordinate axis")

        PROPERTY_DOUBLE(massColumnDensity, "the mass column density along this axis")
        ATTRIBUTE_QUANTITY(massColumnDensity, "masssurfacedensity")
        ATTRIBUTE_MIN_VALUE(massColumnDensity, "]0")

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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NUMBERMATERIALNORMALIZATION_HPP
#define NUMBERMATERIALNORMALIZATION_HPP

#include "MaterialNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A NumberMaterialNormalization object normalizes the amount of material in a geometric medium by
    specifying the total number of entities. */
class NumberMaterialNormalization : public MaterialNormalization
{
    ITEM_CONCRETE(NumberMaterialNormalization, MaterialNormalization,
                  "normalization by defining the total number of entities")

        PROPERTY_DOUBLE(number, "the total number of entities in the material")
        ATTRIBUTE_MIN_VALUE(number, "]0")

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

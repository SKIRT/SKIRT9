/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MASSMATERIALNORMALIZATION_HPP
#define MASSMATERIALNORMALIZATION_HPP

#include "MaterialNormalization.hpp"

////////////////////////////////////////////////////////////////////

/** A MassMaterialNormalization object normalizes the amount of material in a geometric medium by
    specifying the total mass. */
class MassMaterialNormalization : public MaterialNormalization
{
    ITEM_CONCRETE(MassMaterialNormalization, MaterialNormalization, "normalization by defining the total mass")

        PROPERTY_DOUBLE(mass, "the total mass of the material")
        ATTRIBUTE_QUANTITY(mass, "mass")
        ATTRIBUTE_MIN_VALUE(mass, "]0")

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

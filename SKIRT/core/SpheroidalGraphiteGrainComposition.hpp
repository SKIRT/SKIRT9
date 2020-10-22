/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHEROIDALGRAPHITEGRAINCOMPOSITION_HPP
#define SPHEROIDALGRAPHITEGRAINCOMPOSITION_HPP

#include "PolarizedGraphiteGrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The PolarizedGraphiteGrainComposition class represents the optical and calorimetric properties
    of spheroidal silicate dust grains with partial support for polarization. More precisely, the
    current implementation supports polarized thermal emission by (partially) aligned spheroidal
    grains, but assumes spherical grains for scattering and absorption interactions.

    The optical scattering and absorption properties and the calorimetric properties are taken from
    the PolarizedGraphiteGrainComposition class, from which this class derives. The optical
    properties driving the polarization signature for therrmal emission are obtained from
    additional built-in tables or can be provided by the user, as described below.

    TO DO: describe built-in and file input options. */
class SpheroidalGraphiteGrainComposition : public PolarizedGraphiteGrainComposition
{
    ITEM_CONCRETE(SpheroidalGraphiteGrainComposition, PolarizedGraphiteGrainComposition,
                  "a spheroidal graphite dust grain composition with support for polarization")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpheroidalSilicateGrainComposition, "Spheroidal")

        PROPERTY_STRING(spheroidalEmissionTable,
                        "the name of the file tabulating properties for polarized emission by spheroidal grains")

    ITEM_END()

public:
    /** This function returns a brief human-readable identifier for this grain composition. */
    string name() const override;

    /** This function returns information on the resources required for implementing thermal
        emission from aligned spheriodal grains. For more information, see
        GrainComposition::resourcesForSpheroidalEmission. */
    bool resourcesForSpheroidalEmission(bool& resource, double& interpol, string& tableName1,
                                        string& tableName2) const override;
};

////////////////////////////////////////////////////////////////////

#endif

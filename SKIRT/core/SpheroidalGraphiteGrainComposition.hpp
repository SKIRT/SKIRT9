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
    properties driving the polarization signature for thermal emission are obtained from additional
    built-in tables or can be provided by the user, as described below.

    In the current implementation, the internally used optical properties table always assumes a
    fixed spheroidal dust model with a constant alignment. The table generally is a linear
    combination of a table for non-aligned grains and one for aligned grains, combined using a
    linear alignment fraction between 0 and 1. Contrary to the SpheroidalSilicateGrainComposition
    case, the builtin version assumes that graphite grains are non-aligned, and hence consists of a
    single table without an alignment fraction. This builtin table is simply a copy of the
    DraineGraphiteOpticalProps table that contains the additional columns required to handle the
    angular dependence of the optical properties (which is constant in the non-aligned case).

    If this default model is not sufficient, users can provide their own custom tables in SKIRT
    stored table format, e.g. generated using CosTuuM. There are two options: either the user
    computes the tables for a specific alignment fraction and provides a single table, or the user
    provides separate tables for perfectly aligned and non-aligned grains. In the latter case,
    SKIRT will use the alignment fraction to appropriately interpolate between the two tables, as
    in the builtin case.

    The choice between the 3 different scenarios (builtin table without interpolation, a single
    custom table without interpolation or two custom tables with interpolation) is configured
    through an enum. */
class SpheroidalGraphiteGrainComposition : public PolarizedGraphiteGrainComposition
{

    ENUM_DEF(TableType, Builtin, OneTable, TwoTables)
        ENUM_VAL(TableType, Builtin, "builtin resources")
        ENUM_VAL(TableType, OneTable, "single custom table")
        ENUM_VAL(TableType, TwoTables, "two custom tables with interpolation")
    ENUM_END()

    ITEM_CONCRETE(SpheroidalGraphiteGrainComposition, PolarizedGraphiteGrainComposition,
                  "a spheroidal graphite dust grain composition with support for polarization")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpheroidalGraphiteGrainComposition, "Spheroidal")

        PROPERTY_ENUM(tableType, TableType, "the type of emission tables to use")
        ATTRIBUTE_DEFAULT_VALUE(tableType, "Builtin")

        PROPERTY_STRING(emissionTable, "the name of the file tabulating properties for polarized emission by "
                                       "arbitrarily aligned spheroidal grains")
        ATTRIBUTE_RELEVANT_IF(emissionTable, "tableTypeOneTable")

        PROPERTY_STRING(
            alignedEmissionTable,
            "the name of the file tabulating properties for polarized emission by perfectly aligned spheroidal grains")
        ATTRIBUTE_RELEVANT_IF(alignedEmissionTable, "tableTypeTwoTables")

        PROPERTY_STRING(
            nonAlignedEmissionTable,
            "the name of the file tabulating properties for polarized emission by non-aligned spheroidal grains")
        ATTRIBUTE_RELEVANT_IF(nonAlignedEmissionTable, "tableTypeTwoTables")

        PROPERTY_DOUBLE(alignmentFraction,
                        "the alignment fraction of the spheroidal grains with the local magnetic field")
        ATTRIBUTE_DEFAULT_VALUE(alignmentFraction, "1.")
        ATTRIBUTE_MIN_VALUE(alignmentFraction, "0.")
        ATTRIBUTE_MAX_VALUE(alignmentFraction, "1.")
        ATTRIBUTE_RELEVANT_IF(alignmentFraction, "tableTypeTwoTables")

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

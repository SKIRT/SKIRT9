/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHEROIDALSILICATEGRAINCOMPOSITION_HPP
#define SPHEROIDALSILICATEGRAINCOMPOSITION_HPP

#include "PolarizedSilicateGrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The SpheroidalSilicateGrainComposition class represents the optical and calorimetric properties
    of spheroidal silicate dust grains with partial support for polarization. More precisely, the
    current implementation supports polarized thermal emission by (partially) aligned spheroidal
    grains, but assumes spherical grains for scattering and absorption interactions.

    The optical scattering and absorption properties and the calorimetric properties are taken from
    the PolarizedSilicateGrainComposition class, from which this class derives. The optical
    properties driving the polarization signature for thermal emission are obtained from additional
    built-in tables or can be provided by the user, as described below.

    In the current implementation, the internally used optical properties table always assumes a
    fixed spheroidal dust model with a constant alignment. The table generally is a linear
    combination of a table for non-aligned grains and one for aligned grains, combined using a
    linear alignment fraction between 0 and 1. The builtin version of the tables further assumes
    that only grains with sizes larger than 0.1 micron align with the magnetic field, and that the
    spheroidal grains have shapes distributed according to a CDE2 shape distribution (this is the
    fiducial model presented in Vandenbroucke, Baes & Camps, 2020).

    If this fiducial model is not sufficient, users can provide their own custom tables in SKIRT
    stored table format, e.g. generated using CosTuuM. There are two options: either the user
    computes the tables for a specific alignment fraction and provides a single table, or the user
    provides separate tables for perfectly aligned and non-aligned grains. In the latter case,
    SKIRT will use the alignment fraction to appropriately interpolate between the two tables, as
    in the builtin case.

    The choice between the 3 different scenarios (builtin tables with interpolation, a single
    custom table without interpolation or two custom tables with interpolation) is configured
    through an enum. */
class SpheroidalSilicateGrainComposition : public PolarizedSilicateGrainComposition
{

    ENUM_DEF(TableType, Builtin, OneTable, TwoTables)
        ENUM_VAL(TableType, Builtin, "builtin resources")
        ENUM_VAL(TableType, OneTable, "single custom table")
        ENUM_VAL(TableType, TwoTables, "two custom tables with interpolation")
    ENUM_END()

    ITEM_CONCRETE(SpheroidalSilicateGrainComposition, PolarizedSilicateGrainComposition,
                  "a spheroidal silicate dust grain composition with support for polarization")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpheroidalSilicateGrainComposition, "Spheroidal")

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
        ATTRIBUTE_RELEVANT_IF(alignmentFraction, "tableTypeBuiltin|tableTypeTwoTables")

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

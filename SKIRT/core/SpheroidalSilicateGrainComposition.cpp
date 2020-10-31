/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpheroidalSilicateGrainComposition.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

string SpheroidalSilicateGrainComposition::name() const
{
    return "Spheroidal_Polarized_Draine_Silicate";
}

////////////////////////////////////////////////////////////////////

bool SpheroidalSilicateGrainComposition::resourcesForSpheroidalEmission(bool& resource, double& interpol,
                                                                        std::string& tableName1,
                                                                        std::string& tableName2) const
{
    switch (_tableType)
    {
        case TableType::Builtin:
            resource = true;
            interpol = _alignmentFraction;
            tableName1 = "SpheroidalSilicateNonAlignedEmissionOpticalProps";
            tableName2 = "SpheroidalSilicateAlignedEmissionOpticalProps";
            break;
        case TableType::OneTable:
            resource = false;
            interpol = 0;
            tableName1 = _emissionTable;
            tableName2 = string();
            break;
        case TableType::TwoTables:
            resource = false;
            interpol = _alignmentFraction;
            tableName1 = _nonAlignedEmissionTable;
            tableName2 = _alignedEmissionTable;
            break;
    }
    return true;
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpheroidalGraphiteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

string SpheroidalGraphiteGrainComposition::name() const
{
    return "Spheroidal_Polarized_Draine_Graphite";
}

//////////////////////////////////////////////////////////////////////

bool SpheroidalGraphiteGrainComposition::resourcesForSpheroidalEmission(bool& resource, double& interpol,
                                                                        std::string& tableName1,
                                                                        std::string& tableName2) const
{
    resource = false;
    interpol = 0.;
    tableName1 = _spheroidalEmissionTable;
    tableName2 = string();
    return true;
}

//////////////////////////////////////////////////////////////////////

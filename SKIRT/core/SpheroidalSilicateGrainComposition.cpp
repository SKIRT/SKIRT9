/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpheroidalSilicateGrainComposition.hpp"

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
    resource = false;
    interpol = 0.;
    tableName1 = _spheroidalEmissionTable;
    tableName2 = string();
    return true;
}

////////////////////////////////////////////////////////////////////

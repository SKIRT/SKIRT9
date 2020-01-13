/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustEmGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DustEmGrainComposition::DustEmGrainComposition(SimulationItem* parent, GrainType grainType, double bulkMassDensity)
{
    parent->addChild(this);
    _grainType = grainType;
    _bulkMassDensity = bulkMassDensity;
    setup();
}

//////////////////////////////////////////////////////////////////////

namespace
{
    string nameForGrainType(DustEmGrainComposition::GrainType grainType)
    {
        switch (grainType)
        {
            case DustEmGrainComposition::GrainType::aSil: return "aSil";
            case DustEmGrainComposition::GrainType::Gra: return "Gra";
            case DustEmGrainComposition::GrainType::PAH0DL07: return "PAH0_DL07";
            case DustEmGrainComposition::GrainType::PAH1DL07: return "PAH1_DL07";
            case DustEmGrainComposition::GrainType::PAH0MC10: return "PAH0_MC10";
            case DustEmGrainComposition::GrainType::PAH1MC10: return "PAH1_MC10";
            case DustEmGrainComposition::GrainType::CM20: return "CM20";
            case DustEmGrainComposition::GrainType::aOlM5: return "aOlM5";
            case DustEmGrainComposition::GrainType::aPyM5: return "aPyM5";
        }
        return string();
    }
}

//////////////////////////////////////////////////////////////////////

string DustEmGrainComposition::name() const
{
    return "DustEM_" + nameForGrainType(_grainType);
}

//////////////////////////////////////////////////////////////////////

double DustEmGrainComposition::bulkDensity() const
{
    return _bulkMassDensity;
}

//////////////////////////////////////////////////////////////////////

string DustEmGrainComposition::resourceNameForOpticalProps() const
{
    return "DustEM_" + nameForGrainType(_grainType) + "_OpticalProps";
}

//////////////////////////////////////////////////////////////////////

string DustEmGrainComposition::resourceNameForEnthalpies() const
{
    return "DustEM_" + nameForGrainType(_grainType) + "_Enthalpies";
}

//////////////////////////////////////////////////////////////////////

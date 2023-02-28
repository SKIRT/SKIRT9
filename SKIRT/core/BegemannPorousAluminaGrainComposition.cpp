/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BegemannPorousAluminaGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

BegemannPorousAluminaGrainComposition::BegemannPorousAluminaGrainComposition(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

string BegemannPorousAluminaGrainComposition::name() const
{
    return "Begemann_Porous_Alumina";
}

//////////////////////////////////////////////////////////////////////

double BegemannPorousAluminaGrainComposition::bulkDensity() const
{
    return 4.02e3;
}

//////////////////////////////////////////////////////////////////////

string BegemannPorousAluminaGrainComposition::resourceNameForOpticalProps() const
{
    return "BegemannPorousAluminaOpticalProps";
}

//////////////////////////////////////////////////////////////////////

string BegemannPorousAluminaGrainComposition::resourceNameForEnthalpies() const
{
    return "DustEM_aSil_Enthalpies";
}

//////////////////////////////////////////////////////////////////////

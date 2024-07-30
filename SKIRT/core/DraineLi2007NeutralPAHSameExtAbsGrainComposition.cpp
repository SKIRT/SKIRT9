/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineLi2007NeutralPAHSameExtAbsGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DraineLi2007NeutralPAHSameExtAbsGrainComposition::DraineLi2007NeutralPAHSameExtAbsGrainComposition(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007NeutralPAHSameExtAbsGrainComposition::name() const
{
    return "Draine_Neutral_PAH";
}

//////////////////////////////////////////////////////////////////////

double DraineLi2007NeutralPAHSameExtAbsGrainComposition::bulkDensity() const
{
    return 2.24e3;
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007NeutralPAHSameExtAbsGrainComposition::resourceNameForOpticalProps() const
{
    return "DustEM_PAH0_DL07_OpticalProps_SameExtAbs";
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007NeutralPAHSameExtAbsGrainComposition::resourceNameForEnthalpies() const
{
    return "DustEM_PAH0_DL07_Enthalpies";
}

//////////////////////////////////////////////////////////////////////

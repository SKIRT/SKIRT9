/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineLi2007GraphiteSameExtAbsGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DraineLi2007GraphiteSameExtAbsGrainComposition::DraineLi2007GraphiteSameExtAbsGrainComposition(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007GraphiteSameExtAbsGrainComposition::name() const
{
    return "DraineLi2007_Graphite";
}

//////////////////////////////////////////////////////////////////////

double DraineLi2007GraphiteSameExtAbsGrainComposition::bulkDensity() const
{
    return 2.24e3;
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007GraphiteSameExtAbsGrainComposition::resourceNameForOpticalProps() const
{
    return "DustEM_Gra_OpticalProps_SameExtAbs";
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007GraphiteSameExtAbsGrainComposition::resourceNameForEnthalpies() const
{
    return "DustEM_Gra_Enthalpies";
}

//////////////////////////////////////////////////////////////////////

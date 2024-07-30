/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DraineLi2007IonizedPAHSameExtAbsGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

DraineLi2007IonizedPAHSameExtAbsGrainComposition::DraineLi2007IonizedPAHSameExtAbsGrainComposition(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007IonizedPAHSameExtAbsGrainComposition::name() const
{
    return "Draine_Ionized_PAH";
}

//////////////////////////////////////////////////////////////////////

double DraineLi2007IonizedPAHSameExtAbsGrainComposition::bulkDensity() const
{
    return 2.24e3;
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007IonizedPAHSameExtAbsGrainComposition::resourceNameForOpticalProps() const
{
    return "DustEM_PAH1_DL07_OpticalProps_SameExtAbs";
}

//////////////////////////////////////////////////////////////////////

string DraineLi2007IonizedPAHSameExtAbsGrainComposition::resourceNameForEnthalpies() const
{
    return "DustEM_PAH1_DL07_Enthalpies";
}

//////////////////////////////////////////////////////////////////////

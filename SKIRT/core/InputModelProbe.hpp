/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INPUTMODELPROBE_HPP
#define INPUTMODELPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** InputModelProbe is as stub that serves testing purposes. It will be replaced by an actual
    implementation later. */
class InputModelProbe : public Probe
{
    ITEM_CONCRETE(InputModelProbe, Probe, "an input model probe")
        ATTRIBUTE_TYPE_DISPLAYED_IF(InputModelProbe, "Level3")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif

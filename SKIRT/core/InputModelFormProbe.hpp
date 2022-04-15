/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INPUTMODELFORMPROBE_HPP
#define INPUTMODELFORMPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** InputModelFormProbe is as stub that serves testing purposes. It will be replaced by an actual
    implementation later. */
class InputModelFormProbe : public Probe
{
    ITEM_CONCRETE(InputModelFormProbe, Probe, "an input model form probe")
        ATTRIBUTE_TYPE_DISPLAYED_IF(InputModelFormProbe, "Level3")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif

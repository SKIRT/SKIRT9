/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECIALTYPROBE_HPP
#define SPECIALTYPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** SpecialtyProbe is a base class for probes that do \em not cooperate with forms (see the
    ProbeFormBridge class for more information). At this time, this class has no function other
    than separating specialty probes from form probes in the class hierarchy. */
class SpecialtyProbe : public Probe
{
    ITEM_ABSTRACT(SpecialtyProbe, Probe, "a specialty probe")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

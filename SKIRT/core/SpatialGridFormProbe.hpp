/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDFORMPROBE_HPP
#define SPATIALGRIDFORMPROBE_HPP

#include "Form.hpp"
#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** SpatialGridFormProbe is a base class for probes that cooperate with any Form subclass to
    describe how the quantity should be probed, including spatial grid-specific forms and generic
    forms (see the ProbeFormBridge class for more information). The class has no functionality
    other than offering a property for configuring the form. */
class SpatialGridFormProbe : public Probe
{
    ITEM_ABSTRACT(SpatialGridFormProbe, Probe, "a spatial grid form probe")

        PROPERTY_ITEM(form, Form, "the form describing how this quantity should be probed")
        ATTRIBUTE_DEFAULT_VALUE(form, "DefaultCutsForm")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

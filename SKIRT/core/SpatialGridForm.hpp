/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDFORM_HPP
#define SPATIALGRIDFORM_HPP

#include "Form.hpp"

//////////////////////////////////////////////////////////////////////

/** SpatialGridForm is the abstract base class for probe form classes that can be associated only
    with spatial grid probes because they rely on the fact that the probed quantity is discretized
    on the simulation's spatial grid. Refer to the ProbeFormBridge class for more information. */
class SpatialGridForm : public Form
{
    ITEM_ABSTRACT(SpatialGridForm, Form, "a spatial grid-specific probe form")
    ITEM_END()
};

//////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FORM_HPP
#define FORM_HPP

#include "SimulationItem.hpp"
class ProbeFormBridge;

//////////////////////////////////////////////////////////////////////

/** Form is the abstract base class for simulation item classes that describe how a given quantity
    should be probed. Refer to the ProbeFormBridge class for more information. */
class Form : public SimulationItem
{
    ITEM_ABSTRACT(Form, SimulationItem, "a probe form")
    ITEM_END()

public:
    /** This function causes the form to output one or more files for the quantity being probed
        according to the information provided by the specified ProbeFormBridge instance. The
        function must be implemented by each subclass. */
    virtual void writeQuantity(const ProbeFormBridge* bridge) const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif

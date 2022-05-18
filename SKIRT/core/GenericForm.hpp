/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GENERICFORM_HPP
#define GENERICFORM_HPP

#include "Form.hpp"

//////////////////////////////////////////////////////////////////////

/** GenericForm is the abstract base class for probe form classes that can be associated with any
    probe, including both spatial grid and input model probes. Refer to the ProbeFormBridge class
    for more information. */
class GenericForm : public Form
{
    ITEM_ABSTRACT(GenericForm, Form, "a generic probe form")
    ITEM_END()
};

//////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProbeFormBridge.hpp"

////////////////////////////////////////////////////////////////////

ProbeFormBridge::ProbeFormBridge(const Probe* probe, const Form* form)
{
    _probe = probe;
    _form = form;
}

////////////////////////////////////////////////////////////////////

void ProbeFormBridge::setQuantity(int numValues, string quantity, string projectedQuantity, string description,
                                  string projectedDescription, AddColumnHeaders addColumnHeaders)
{
    _numValues = numValues;
    _quantity = quantity;
    _projectedQuantity = projectedQuantity;
    _description = description;
    _projectedDescription = projectedDescription;
    _addColumnHeaders = addColumnHeaders;
}

////////////////////////////////////////////////////////////////////

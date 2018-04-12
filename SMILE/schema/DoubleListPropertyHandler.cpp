/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DoubleListPropertyHandler.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"

////////////////////////////////////////////////////////////////////

bool DoubleListPropertyHandler::hasDefaultValue() const
{
    return isValidDoubleList(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

void DoubleListPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

vector<double> DoubleListPropertyHandler::defaultValue() const
{
    return toDoubleList(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

vector<double> DoubleListPropertyHandler::value() const
{
    return target()->getDoubleListProperty(property());
}

////////////////////////////////////////////////////////////////////

void DoubleListPropertyHandler::setValue(vector<double> value)
{
    target()->setDoubleListProperty(property(), value);
    setChanged();
}

////////////////////////////////////////////////////////////////////

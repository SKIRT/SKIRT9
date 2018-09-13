/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DoubleListPropertyHandler.hpp"
#include "Item.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"

////////////////////////////////////////////////////////////////////

bool DoubleListPropertyHandler::isValidValue(string value) const
{
    return isValidDoubleList(value);
}

////////////////////////////////////////////////////////////////////

void DoubleListPropertyHandler::insertNames()
{
    if (!value().empty()) nameManager()->insert(property()->name());
    nameManager()->insertFromConditionalValue(property()->insert());
}

////////////////////////////////////////////////////////////////////

void DoubleListPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

vector<double> DoubleListPropertyHandler::defaultValue() const
{
    return toDoubleList(nameManager()->evaluateConditionalValue(property()->defaultValue()));
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

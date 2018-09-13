/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StringPropertyHandler.hpp"
#include "Item.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool StringPropertyHandler::isValidValue(string value) const
{
    return !value.empty();
}

////////////////////////////////////////////////////////////////////

void StringPropertyHandler::insertNames()
{
    if (!value().empty()) nameManager()->insert(property()->name());
    nameManager()->insertFromConditionalValue(property()->insert());
}

////////////////////////////////////////////////////////////////////

void StringPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

string StringPropertyHandler::defaultValue() const
{
    return nameManager()->evaluateConditionalValue(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

string StringPropertyHandler::value() const
{
    return target()->getStringProperty(property());
}

////////////////////////////////////////////////////////////////////

void StringPropertyHandler::setValue(string value)
{
    target()->setStringProperty(property(), value);
    setChanged();
}

////////////////////////////////////////////////////////////////////

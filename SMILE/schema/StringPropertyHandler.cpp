/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StringPropertyHandler.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool StringPropertyHandler::isTrueInCondition() const
{
    return !value().empty();
}

////////////////////////////////////////////////////////////////////

bool StringPropertyHandler::isOptional() const
{
    return StringUtils::toBool(property()->optional()) && property()->defaultValue().empty();
}

////////////////////////////////////////////////////////////////////

bool StringPropertyHandler::hasDefaultValue() const
{
    return !property()->defaultValue().empty();
}

////////////////////////////////////////////////////////////////////

void StringPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

string StringPropertyHandler::defaultValue() const
{
    return property()->defaultValue();
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

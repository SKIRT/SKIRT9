/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoolPropertyHandler.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool BoolPropertyHandler::isTrueInCondition() const
{
    return value();
}

////////////////////////////////////////////////////////////////////

bool BoolPropertyHandler::hasDefaultValue() const
{
    return StringUtils::isValidBool(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

void BoolPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

bool BoolPropertyHandler::defaultValue() const
{
    return StringUtils::toBool(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

bool BoolPropertyHandler::value() const
{
    return target()->getBoolProperty(property());
}

////////////////////////////////////////////////////////////////////

void BoolPropertyHandler::setValue(bool value)
{
    target()->setBoolProperty(property(), value);
    setChanged();
}

////////////////////////////////////////////////////////////////////

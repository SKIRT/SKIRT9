/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoolPropertyHandler.hpp"
#include "Item.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool BoolPropertyHandler::isValidValue(string value) const
{
    return StringUtils::isValidBool(value);
}

////////////////////////////////////////////////////////////////////

void BoolPropertyHandler::insertNames()
{
    if (value()) nameManager()->insert(property()->name());
    nameManager()->insertFromConditionalValue(property()->insert());
}

////////////////////////////////////////////////////////////////////

void BoolPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

bool BoolPropertyHandler::defaultValue() const
{
    return StringUtils::toBool(nameManager()->evaluateConditionalValue(property()->defaultValue()));
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

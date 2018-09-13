/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EnumPropertyHandler.hpp"
#include "Item.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool EnumPropertyHandler::isValidValue(string value) const
{
    return StringUtils::contains(property()->enumNames(), value);
}

////////////////////////////////////////////////////////////////////

void EnumPropertyHandler::insertNames()
{
    nameManager()->insert(property()->name() + value());
    nameManager()->insertFromConditionalValue(property()->insert());
}

////////////////////////////////////////////////////////////////////

void EnumPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

vector<string> EnumPropertyHandler::values() const
{
    return property()->enumNames();
}

////////////////////////////////////////////////////////////////////

vector<string> EnumPropertyHandler::titlesForValues() const
{
    return property()->enumTitles();
}

////////////////////////////////////////////////////////////////////

string EnumPropertyHandler::defaultValue() const
{
    return nameManager()->evaluateConditionalValue(property()->defaultValue());
}

////////////////////////////////////////////////////////////////////

string EnumPropertyHandler::value() const
{
    return target()->getEnumProperty(property());
}

////////////////////////////////////////////////////////////////////

string EnumPropertyHandler::titleForValue() const
{
    return property()->enumTitle(value());
}

////////////////////////////////////////////////////////////////////

void EnumPropertyHandler::setValue(string value)
{
    target()->setEnumProperty(property(), value);
    setChanged();
}

////////////////////////////////////////////////////////////////////

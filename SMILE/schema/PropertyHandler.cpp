/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PropertyHandler.hpp"
#include "BooleanExpression.hpp"
#include "Item.hpp"
#include "ItemUtils.hpp"
#include "NameManager.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

PropertyHandler::PropertyHandler(Item* target, const PropertyDef* property,
                                 const SchemaDef* schema, NameManager* nameMgr)
    : _target(target), _property(property), _schema(schema), _nameMgr(nameMgr)
{
}

////////////////////////////////////////////////////////////////////

string PropertyHandler::type() const
{
    return target()->type();
}

////////////////////////////////////////////////////////////////////

string PropertyHandler::name() const
{
    return property()->name();
}

////////////////////////////////////////////////////////////////////

string PropertyHandler::title() const
{
    return property()->title();
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isRelevant() const
{
    return nameManager()->evaluateBoolean(property()->relevantIf());
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isDisplayed() const
{
    return nameManager()->evaluateBoolean(property()->displayedIf());
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isRequired() const
{
    return nameManager()->evaluateBoolean(property()->requiredIf());
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::hasDefaultValue() const
{
    return isValidValue(nameManager()->evaluateConditionalValue(property()->defaultValue()));
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isCompound() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::hasChanged() const
{
    return _changed;
}

////////////////////////////////////////////////////////////////////

void PropertyHandler::setChanged()
{
    _changed = true;
}

////////////////////////////////////////////////////////////////////

void PropertyHandler::setConfigured(bool configured)
{
    ItemUtils::setPropertyConfigured(target(), name(), configured);
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isConfigured()
{
    return ItemUtils::isPropertyConfigured(target(), name());
}

////////////////////////////////////////////////////////////////////

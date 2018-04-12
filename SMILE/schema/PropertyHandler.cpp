/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PropertyHandler.hpp"
#include "BooleanExpression.hpp"
#include "Item.hpp"
#include "ItemUtils.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

PropertyHandler::PropertyHandler(Item* target, const PropertyDef* property, const SchemaDef* schema)
    : _target(target), _property(property), _schema(schema)
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

bool PropertyHandler::isSilent() const
{
    return StringUtils::toBool(property()->silent());
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::hasRelevantIf() const
{
    return !property()->relevantIf().empty();
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isRelevant() const
{
    // get the expression
    string expression = property()->relevantIf();
    if (expression.empty()) return true;

    // evaluate the expression
    return BooleanExpression::evaluate(expression, [this] (string name)
        {
            // construct a handler for the target property and evaluate our relevancy
            auto handler = schema()->createPropertyHandler(target(), name);
            return handler->isRelevant() && handler->isTrueInCondition();
        });
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isTrueInCondition() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isOptional() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::hasDefaultValue() const
{
    return false;
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

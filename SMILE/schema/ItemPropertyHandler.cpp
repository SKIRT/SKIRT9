/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ItemPropertyHandler.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

bool ItemPropertyHandler::isTrueInCondition() const
{
    return value() != 0;
}

////////////////////////////////////////////////////////////////////

void ItemPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

Item* ItemPropertyHandler::value() const
{
    return target()->getItemProperty(property());
}

////////////////////////////////////////////////////////////////////

bool ItemPropertyHandler::setValue(Item* value)
{
    if (schema()->inherits(value->type(), baseType()))
    {
        target()->setItemProperty(property(), value);
        setChanged();
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

bool ItemPropertyHandler::setToNewItemOfType(string type)
{
    if (schema()->inherits(type, baseType()))
    {
        target()->setItemProperty(property(), schema()->createItem(type).release());
        setChanged();
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

void ItemPropertyHandler::setToNull()
{
    target()->setItemProperty(property(), nullptr);
    setChanged();
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ItemPropertyHandler.hpp"
#include "Item.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

vector<Item*> ItemPropertyHandler::children() const
{
    vector<Item*> result;
    if (value()) result.push_back(value());
    return result;
}

////////////////////////////////////////////////////////////////////

void ItemPropertyHandler::insertNames()
{
    if (value()) nameManager()->insert(property()->name());
    nameManager()->insertFromConditionalValue(property()->insert());

    if (value())
    {
        nameManager()->insert(schema()->ascendants(value()->type()));
        nameManager()->insertFromConditionalValue(schema()->toBeInserted(value()->type()));
    }
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
    if (isValidValue(value->type()))
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
    if (isValidValue(type))
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

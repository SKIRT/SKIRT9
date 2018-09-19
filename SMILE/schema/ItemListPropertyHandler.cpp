/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ItemListPropertyHandler.hpp"
#include "Item.hpp"
#include "ItemUtils.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

vector<Item*> ItemListPropertyHandler::children() const
{
    return value();
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyHandler::insertNames()
{
    if (!value().empty()) nameManager()->insert(property()->name());
    nameManager()->insertFromConditionalValue(property()->insert());

    for (auto item : value())
    {
        nameManager()->insert(schema()->ascendants(item->type()));
        nameManager()->insertFromConditionalValue(schema()->toBeInserted(item->type()));
    }
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyHandler::acceptVisitor(PropertyHandlerVisitor* visitor)
{
    visitor->visitPropertyHandler(this);
}

////////////////////////////////////////////////////////////////////

vector<Item*> ItemListPropertyHandler::value() const
{
    return target()->getItemListProperty(property());
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyHandler::setToEmpty()
{
    target()->clearItemListProperty(property());
}

////////////////////////////////////////////////////////////////////

bool ItemListPropertyHandler::addValue(Item* value)
{
    if (isValidValue(value->type()))
    {
        target()->insertIntoItemListProperty(property(), -1, value);
        setChanged();
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

bool ItemListPropertyHandler::addNewItemOfType(string type)
{
    if (isValidValue(type))
    {
        target()->insertIntoItemListProperty(property(), -1, schema()->createItem(type).release());
        setChanged();
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

bool ItemListPropertyHandler::insertValue(int index, Item* value)
{
    if (isValidValue(value->type()))
    {
        target()->insertIntoItemListProperty(property(), index, value);
        setChanged();
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

bool ItemListPropertyHandler::insertNewItemOfType(int index, string type)
{
    if (isValidValue(type))
    {
        target()->insertIntoItemListProperty(property(), index, schema()->createItem(type).release());
        setChanged();
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

bool ItemListPropertyHandler::removeValueAt(int index)
{
    if (index >= 0 && static_cast<size_t>(index) < value().size())
    {
        target()->removeFromItemListProperty(property(), index);
        setChanged();
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

void ItemListPropertyHandler::storeSelectedRow(int row)
{
    ItemUtils::storeSelectedRow(target(), name(), row);
}

////////////////////////////////////////////////////////////////////

int ItemListPropertyHandler::retrieveSelectedRow()
{
    return ItemUtils::retrieveSelectedRow(target(), name());
}

////////////////////////////////////////////////////////////////////

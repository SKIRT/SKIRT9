/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ItemUtils.hpp"
#include "FatalError.hpp"
#include "Item.hpp"
#include "SchemaDef.hpp"

////////////////////////////////////////////////////////////////////

void ItemUtils::setPropertyConfiguredState(Item* item, string property, int configured)
{
    item->setUtilityProperty(property + "@configured", configured);
}

////////////////////////////////////////////////////////////////////

void ItemUtils::setHierarchyConfigured(const SchemaDef* schema, Item* root)
{
    // process all immediate properties of the specified root item
    for (auto property : schema->properties(root->type())) setPropertyConfiguredState(root, property, 1);

    // process all children of the specified root item
    for (auto child : root->children()) setHierarchyConfigured(schema, child);
}

////////////////////////////////////////////////////////////////////

int ItemUtils::propertyConfiguredState(Item* item, string property)
{
    try
    {
        return item->getUtilityProperty(property + "@configured");
    }
    catch (const FatalError&)
    {}
    return 0;
}

////////////////////////////////////////////////////////////////////

void ItemUtils::setItemComplete(Item* item)
{
    item->setUtilityProperty("item@complete", true);
}

////////////////////////////////////////////////////////////////////

void ItemUtils::setHierarchyComplete(Item* root)
{
    // process the specified root item
    setItemComplete(root);

    // process all children of the specified root item
    for (auto child : root->children()) setHierarchyComplete(child);
}

////////////////////////////////////////////////////////////////////

void ItemUtils::setItemIncomplete(Item* item)
{
    while (item)
    {
        item->setUtilityProperty("item@complete", false);
        item = item->parent();
    }
}

////////////////////////////////////////////////////////////////////

bool ItemUtils::isItemComplete(Item* item)
{
    try
    {
        return item->getUtilityProperty("item@complete") ? true : false;
    }
    catch (const FatalError&)
    {}
    return false;
}

////////////////////////////////////////////////////////////////////

void ItemUtils::storeSelectedRow(Item* item, string property, int row)
{
    item->setUtilityProperty(property + "@row", row);
}

////////////////////////////////////////////////////////////////////

int ItemUtils::retrieveSelectedRow(Item* item, string property)
{
    try
    {
        return item->getUtilityProperty(property + "@row");
    }
    catch (const FatalError&)
    {}
    return 0;
}

////////////////////////////////////////////////////////////////////

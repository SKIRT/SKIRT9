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

PropertyHandler::PropertyHandler(Item* target, const PropertyDef* property, const SchemaDef* schema,
                                 NameManager* nameMgr)
    : _target(target), _property(property), _schema(schema), _nameMgr(nameMgr)
{}

////////////////////////////////////////////////////////////////////

Item* PropertyHandler::root() const
{
    Item* result = target();
    while (result->parent()) result = result->parent();
    return result;
}

////////////////////////////////////////////////////////////////////

vector<Item*> PropertyHandler::children() const
{
    return vector<Item*>();
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

bool PropertyHandler::isSilent() const
{
    return !isRelevant() || (!isDisplayed() && (!isRequired() || hasDefaultValue()));
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
    insertNames();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function is used by rebuildNames() to recursively insert the names
    // starting at the specified item and ending just before the specified target item/property
    bool insertNamesRecursively(Item* item, const SchemaDef* schema, NameManager* nameMgr, Item* targetItem,
                                string targetProperty)
    {
        nameMgr->pushLocal();

        // loop over all properties of this item
        for (const string& property : schema->properties(item->type()))
        {
            // abort the recursion just before the target property in the target item
            if (item == targetItem && property == targetProperty) return false;

            // insert the names for this property
            auto handler = schema->createPropertyHandler(item, property, nameMgr);
            handler->insertNames();

            // recursively handle the properties of the item's children, if any
            for (Item* child : handler->children())
            {
                if (!insertNamesRecursively(child, schema, nameMgr, targetItem, targetProperty)) return false;
            }
        }

        nameMgr->popLocal();
        return true;
    }
}

////////////////////////////////////////////////////////////////////

void PropertyHandler::rebuildNames()
{
    // clear the name sets
    nameManager()->clearAll();

    // start the recursion with the root item
    insertNamesRecursively(root(), schema(), nameManager(), target(), property()->name());
}

////////////////////////////////////////////////////////////////////

void PropertyHandler::setNotConfigured()
{
    ItemUtils::setPropertyConfiguredState(target(), name(), 0);
}

////////////////////////////////////////////////////////////////////

void PropertyHandler::setConfiguredToDefault()
{
    ItemUtils::setPropertyConfiguredState(target(), name(), 2);
}

////////////////////////////////////////////////////////////////////

void PropertyHandler::setConfiguredByUser(bool valid)
{
    ItemUtils::setPropertyConfiguredState(target(), name(), valid ? 1 : -1);
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isConfiguredByUser()
{
    return abs(ItemUtils::propertyConfiguredState(target(), name())) == 1;
}

////////////////////////////////////////////////////////////////////

bool PropertyHandler::isConfigured()
{
    return ItemUtils::propertyConfiguredState(target(), name()) > 0;
}

////////////////////////////////////////////////////////////////////

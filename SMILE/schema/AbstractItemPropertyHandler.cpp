/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AbstractItemPropertyHandler.hpp"
#include "Item.hpp"
#include "PropertyDef.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"
#include <unordered_set>

////////////////////////////////////////////////////////////////////

bool AbstractItemPropertyHandler::isOptional() const
{
    return StringUtils::toBool(property()->optional());
}

////////////////////////////////////////////////////////////////////

bool AbstractItemPropertyHandler::hasDefaultValue() const
{
    string type = property()->defaultValue();
    return !type.empty() && schema()->inherits(type, baseType());
}

////////////////////////////////////////////////////////////////////

bool AbstractItemPropertyHandler::isCompound() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

string AbstractItemPropertyHandler::baseType() const
{
    return property()->base();
}

////////////////////////////////////////////////////////////////////

string AbstractItemPropertyHandler::defaultType() const
{
    return property()->defaultValue();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // Returns the root of the hierarchy in which the specified item resides.
    Item* getRoot(Item* item)
    {
        while (item->parent()) item = item->parent();
        return item;
    }

    // Adds to the specified set:
    //  - the types of all items present in the specified hierarchy, and
    //  - the types of all their compile-time ascendants.
    // The function calls itself recursively to process the children of the specified root item.
    void addHierarchyTypeNames(const SchemaDef* schema, Item* root, std::unordered_set<string>& keys)
    {
        // process the specified root item
        for (auto key : schema->ascendants(root->type())) keys.insert(key);

        // process all children of the specified root item
        for (auto child : root->children()) addHierarchyTypeNames(schema, child, keys);
    }
}

////////////////////////////////////////////////////////////////////

vector<string> AbstractItemPropertyHandler::allowedDescendants()
{
    std::unordered_set<string> keys;
    addHierarchyTypeNames(schema(), getRoot(target()), keys);
    return schema()->allowedDescendants(baseType(), keys);
}

////////////////////////////////////////////////////////////////////

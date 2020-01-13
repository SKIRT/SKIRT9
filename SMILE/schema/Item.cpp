/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Item.hpp"
#include "FatalError.hpp"
#include "ItemRegistry.hpp"
#include "PropertyAccessor.hpp"
#include "PropertyDef.hpp"
#include "StringUtils.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

// use a custom class to avoid including std::unordered_map in the Item.hpp header
class Item::UtilityData : public std::unordered_map<string, int>
{};

////////////////////////////////////////////////////////////////////

const char* Item::ii_type() const
{
    return "Item";
}

////////////////////////////////////////////////////////////////////

void Item::ii_loadItemInfo()
{
    ItemRegistry::beginType("Item", "", "a SMILE item");
}

////////////////////////////////////////////////////////////////////

Item::~Item()
{
    delete _utility;
    while (_children.size()) delete _children.back();  // child destructor calls removeChild()
    if (_parent) _parent->removeChild(this);
}

////////////////////////////////////////////////////////////////////

void Item::setParent(Item* parent)
{
    _parent = parent;
}

////////////////////////////////////////////////////////////////////

void Item::removeChild(Item* child)
{
    _children.erase(std::find(_children.begin(), _children.end(), child));
}

////////////////////////////////////////////////////////////////////

void Item::addChild(Item* child)
{
    if (child)
    {
        _children.push_back(child);
        child->setParent(this);
    }
}

////////////////////////////////////////////////////////////////////

void Item::destroyChild(Item* child)
{
    delete child;  // child destructor calls removeChild()
}

////////////////////////////////////////////////////////////////////

Item* Item::parent() const
{
    return _parent;
}

////////////////////////////////////////////////////////////////////

const vector<Item*>& Item::children() const
{
    return _children;
}

////////////////////////////////////////////////////////////////////

string Item::type() const
{
    return ii_type();
}

////////////////////////////////////////////////////////////////////

void Item::setStringProperty(const PropertyDef* property, string value)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<string>*>(property->accessor());
    accessor->setTargetValue(this, value);
}

////////////////////////////////////////////////////////////////////

void Item::setBoolProperty(const PropertyDef* property, bool value)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<bool>*>(property->accessor());
    accessor->setTargetValue(this, value);
}

////////////////////////////////////////////////////////////////////

void Item::setIntProperty(const PropertyDef* property, int value)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<int>*>(property->accessor());
    accessor->setTargetValue(this, value);
}

////////////////////////////////////////////////////////////////////

void Item::setEnumProperty(const PropertyDef* property, string value)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<int>*>(property->accessor());
    accessor->setTargetValue(this, StringUtils::indexOf(property->enumNames(), value));
}

////////////////////////////////////////////////////////////////////

void Item::setDoubleProperty(const PropertyDef* property, double value)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<double>*>(property->accessor());
    accessor->setTargetValue(this, value);
}

////////////////////////////////////////////////////////////////////

void Item::setDoubleListProperty(const PropertyDef* property, vector<double> value)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<vector<double>>*>(property->accessor());
    accessor->setTargetValue(this, value);
}

////////////////////////////////////////////////////////////////////

void Item::setItemProperty(const PropertyDef* property, Item* item)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<Item*>*>(property->accessor());
    accessor->setTargetValue(this, item);
}

////////////////////////////////////////////////////////////////////

void Item::clearItemListProperty(const PropertyDef* property)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<vector<Item*>>*>(property->accessor());
    accessor->clearTargetValue(this);
}

////////////////////////////////////////////////////////////////////

void Item::insertIntoItemListProperty(const PropertyDef* property, int index, Item* item)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<vector<Item*>>*>(property->accessor());
    accessor->insertIntoTargetValue(this, index, item);
}

////////////////////////////////////////////////////////////////////

void Item::removeFromItemListProperty(const PropertyDef* property, int index)
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<vector<Item*>>*>(property->accessor());
    accessor->removeFromTargetValue(this, index);
}

////////////////////////////////////////////////////////////////////

string Item::getStringProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<string>*>(property->accessor());
    return accessor->targetValue(this);
}

////////////////////////////////////////////////////////////////////

bool Item::getBoolProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<bool>*>(property->accessor());
    return accessor->targetValue(this);
}

////////////////////////////////////////////////////////////////////

int Item::getIntProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<int>*>(property->accessor());
    return accessor->targetValue(this);
}

////////////////////////////////////////////////////////////////////

string Item::getEnumProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<int>*>(property->accessor());
    return property->enumNames()[accessor->targetValue(this)];
}

////////////////////////////////////////////////////////////////////

double Item::getDoubleProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<double>*>(property->accessor());
    return accessor->targetValue(this);
}

////////////////////////////////////////////////////////////////////

vector<double> Item::getDoubleListProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<vector<double>>*>(property->accessor());
    return accessor->targetValue(this);
}

////////////////////////////////////////////////////////////////////

Item* Item::getItemProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<Item*>*>(property->accessor());
    return accessor->targetValue(this);
}

////////////////////////////////////////////////////////////////////

vector<Item*> Item::getItemListProperty(const PropertyDef* property) const
{
    auto accessor = dynamic_cast<const TypedPropertyAccessor<vector<Item*>>*>(property->accessor());
    return accessor->targetValue(this);
}

////////////////////////////////////////////////////////////////////

void Item::setUtilityProperty(string name, int value)
{
    if (!_utility) _utility = new Item::UtilityData;
    (*_utility)[name] = value;
}

////////////////////////////////////////////////////////////////////

int Item::getUtilityProperty(string name) const
{
    if (_utility && _utility->count(name))
        return _utility->at(name);
    else
        throw FATALERROR("Unknow ghost property " + name);
}

////////////////////////////////////////////////////////////////////

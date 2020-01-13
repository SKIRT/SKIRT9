/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GhostItem.hpp"
#include "FatalError.hpp"
#include "PropertyDef.hpp"

////////////////////////////////////////////////////////////////////

void GhostItem::setStringProperty(const PropertyDef* property, string value)
{
    _stringProperties[property->name()] = value;
}

////////////////////////////////////////////////////////////////////

void GhostItem::setBoolProperty(const PropertyDef* property, bool value)
{
    _boolProperties[property->name()] = value;
}

////////////////////////////////////////////////////////////////////

void GhostItem::setIntProperty(const PropertyDef* property, int value)
{
    _intProperties[property->name()] = value;
}

////////////////////////////////////////////////////////////////////

void GhostItem::setEnumProperty(const PropertyDef* property, string value)
{
    _enumProperties[property->name()] = value;
}

////////////////////////////////////////////////////////////////////

void GhostItem::setDoubleProperty(const PropertyDef* property, double value)
{
    _doubleProperties[property->name()] = value;
}

////////////////////////////////////////////////////////////////////

void GhostItem::setDoubleListProperty(const PropertyDef* property, vector<double> value)
{
    _doubleListProperties[property->name()] = value;
}

////////////////////////////////////////////////////////////////////

void GhostItem::setItemProperty(const PropertyDef* property, Item* item)
{
    auto& prop = _itemProperties[property->name()];
    destroyChild(prop);
    prop = item;
    addChild(prop);
}

////////////////////////////////////////////////////////////////////

void GhostItem::clearItemListProperty(const PropertyDef* property)
{
    auto& list = _itemListProperties[property->name()];
    for (auto item : list) destroyChild(item);
    list.clear();
}

////////////////////////////////////////////////////////////////////

void GhostItem::insertIntoItemListProperty(const PropertyDef* property, int index, Item* item)
{
    auto& list = _itemListProperties[property->name()];
    if (index >= 0 && static_cast<size_t>(index) < list.size())
    {
        list.insert(list.cbegin() + index, item);
    }
    else
    {
        list.push_back(item);
    }
    addChild(item);
}

////////////////////////////////////////////////////////////////////

void GhostItem::removeFromItemListProperty(const PropertyDef* property, int index)
{
    auto& list = _itemListProperties[property->name()];
    if (index >= 0 && static_cast<size_t>(index) < list.size())
    {
        destroyChild(list[index]);
        list.erase(list.cbegin() + index, list.cbegin() + index + 1);
    }
}

////////////////////////////////////////////////////////////////////

string GhostItem::type() const
{
    return _type;
}

////////////////////////////////////////////////////////////////////

string GhostItem::getStringProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_stringProperties.count(name))
        return _stringProperties.at(name);
    else
        throw FATALERROR("Unknow string property " + name);
}

////////////////////////////////////////////////////////////////////

bool GhostItem::getBoolProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_boolProperties.count(name))
        return _boolProperties.at(name);
    else
        throw FATALERROR("Unknow Boolean property " + name);
}

////////////////////////////////////////////////////////////////////

int GhostItem::getIntProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_intProperties.count(name))
        return _intProperties.at(name);
    else
        throw FATALERROR("Unknow integer property " + name);
}

////////////////////////////////////////////////////////////////////

string GhostItem::getEnumProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_enumProperties.count(name))
        return _enumProperties.at(name);
    else
        throw FATALERROR("Unknow enumeration property " + name);
}

////////////////////////////////////////////////////////////////////

double GhostItem::getDoubleProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_doubleProperties.count(name))
        return _doubleProperties.at(name);
    else
        throw FATALERROR("Unknow double property " + name);
}

////////////////////////////////////////////////////////////////////

vector<double> GhostItem::getDoubleListProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_doubleListProperties.count(name))
        return _doubleListProperties.at(name);
    else
        throw FATALERROR("Unknow double list property " + name);
}

////////////////////////////////////////////////////////////////////

Item* GhostItem::getItemProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_itemProperties.count(name))
        return _itemProperties.at(name);
    else
        throw FATALERROR("Unknow item property " + name);
}

////////////////////////////////////////////////////////////////////

vector<Item*> GhostItem::getItemListProperty(const PropertyDef* property) const
{
    string name = property->name();
    if (_itemListProperties.count(name))
        return _itemListProperties.at(name);
    else
        throw FATALERROR("Unknow item list property " + name);
}

////////////////////////////////////////////////////////////////////

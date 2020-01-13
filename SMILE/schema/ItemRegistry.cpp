/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ItemRegistry.hpp"
#include "FatalError.hpp"
#include "Item.hpp"
#include "PropertyAccessor.hpp"
#include "PropertyDef.hpp"
#include "SchemaDef.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

namespace
{
    std::unordered_map<string, SchemaDef> _schemas;
    SchemaDef* _targetSchema{nullptr};
    TypeDef* _targetType{nullptr};
    PropertyDef* _targetProperty{nullptr};
    int _enumCount{-1};
    int _enumIndex{-1};
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::beginSchema(string name, string title, string version, string extension, string root, string type,
                               string format, string url)
{
    // verify that the schema name is not in use
    if (_schemas.count(name)) throw FATALERROR("Schema '" + name + "' is already defined in the registry");

    // create a schema with the specified basic attributes
    _schemas.emplace(std::piecewise_construct, std::forward_as_tuple(name),
                     std::forward_as_tuple(name, title, version, extension, root, type, format, url));

    // remember the new schema as the current target
    _targetSchema = &_schemas.at(name);
    _targetType = nullptr;
    _targetProperty = nullptr;

    // always add the Item base class
    ItemRegistry::add<Item>();
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::addUnitDefImpl(const UnitDef& unitDef)
{
    if (!_targetSchema) throw FATALERROR("Adding unit definition without target schema");
    _targetSchema->loadUnitDef(unitDef);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::finalize()
{
    _targetSchema = nullptr;
    _targetType = nullptr;
    _targetProperty = nullptr;
    _schemas.clear();
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::beginType(const char* type, const char* baseType, const char* title, Instantiator instantiator)
{
    if (!_targetSchema) throw FATALERROR("Adding type without target schema");
    _targetType = &_targetSchema->addTypeDef(type, baseType, title, instantiator);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setTypeAllowedIf(const char* expression)
{
    _targetType->setAllowedIf(expression);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setTypeDisplayedIf(const char* expression)
{
    _targetType->setDisplayedIf(expression);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setTypeInsert(const char* expression)
{
    _targetType->setInsert(expression);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setSubPropertyIndexHere()
{
    _targetType->setSubPropertyIndexHere();
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::beginProperty(const char* type, const char* name, const char* title,
                                 const PropertyAccessor* accessor)
{
    if (!_targetType) throw FATALERROR("Adding property without target type");
    _targetProperty = &_targetType->addPropertyDef(type, name, title);
    _targetProperty->setAccessor(accessor);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::beginEnum(int enumCount)
{
    if (enumCount <= 0) throw FATALERROR("Enumeration type has no enumeration elements");
    _enumIndex = 0;
    _enumCount = enumCount;
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::addEnum(int enumIndex, const char* enumName, const char* enumTitle)
{
    if (enumIndex != _enumIndex) throw FATALERROR("Enumeration element is added out of order");
    _targetProperty->addEnumeration(enumName, enumTitle);
    _enumIndex++;
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::endEnum()
{
    if (_enumIndex != _enumCount) throw FATALERROR("Number of enumeration elements added does not match");
    _enumIndex = -1;
    _enumCount = -1;
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setRelevantIf(const char* expression)
{
    _targetProperty->setRelevantIf(expression);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setDisplayedIf(const char* expression)
{
    _targetProperty->setDisplayedIf(expression);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setRequiredIf(const char* expression)
{
    _targetProperty->setRequiredIf(expression);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setInsert(const char* expression)
{
    _targetProperty->setInsert(expression);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setDefaultValue(const char* value)
{
    _targetProperty->setDefaultValue(value);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setMinValue(const char* value)
{
    _targetProperty->setMinValue(value);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setMaxValue(const char* value)
{
    _targetProperty->setMaxValue(value);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setQuantity(const char* value)
{
    _targetProperty->setQuantity(value);
}

////////////////////////////////////////////////////////////////////

void ItemRegistry::setBase(const char* type)
{
    _targetProperty->setBase(type);
}

////////////////////////////////////////////////////////////////////

const SchemaDef* ItemRegistry::getSchemaDef(string name)
{
    if (_schemas.count(name))
        return &_schemas.at(name);
    else
        throw FATALERROR("Unknow schema definition " + name);
}

////////////////////////////////////////////////////////////////////

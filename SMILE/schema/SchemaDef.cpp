/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SchemaDef.hpp"
#include "BoolPropertyHandler.hpp"
#include "BooleanExpression.hpp"
#include "DoubleListPropertyHandler.hpp"
#include "DoublePropertyHandler.hpp"
#include "EnumPropertyHandler.hpp"
#include "FatalError.hpp"
#include "GhostItem.hpp"
#include "IntPropertyHandler.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemPropertyHandler.hpp"
#include "PropertyDef.hpp"
#include "StringPropertyHandler.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "XmlReader.hpp"
#include "XmlWriter.hpp"

////////////////////////////////////////////////////////////////////

string SchemaDef::getSchemaTitle(string filePath)
{
    try
    {
        // open the file for XML parsing
        XmlReader reader(filePath);

        // read and verify the first and second elements
        if (reader.readNextStartElement() && reader.elementName() == "smile-schema" && reader.readNextStartElement()
            && reader.elementName() == "Schema")
        {
            // get the value of the title attribute (or empty string if there is no such attribute)
            return reader.attributeValue("title");
        }
    }

    // catch XML-reader errors and return an empty string
    catch (const FatalError&)
    {}
    return string();
}

////////////////////////////////////////////////////////////////////

bool SchemaDef::isCompatible(string schemaFilePath, string dataFilePath)
{
    try
    {
        // fields to verify between schema and dataset
        string root1, root2;      // the name of the root element in datasets described by this schema
        string type1, type2;      // the type of the top-level item in datasets described by this schema
        string format1, format2;  // the version of the described data format (specified on the root element)

        // get the value of these fields from the schema file
        {
            XmlReader reader(schemaFilePath);
            if (reader.readNextStartElement() && reader.elementName() == "smile-schema" && reader.readNextStartElement()
                && reader.elementName() == "Schema")
            {
                root1 = reader.attributeValue("root");
                type1 = reader.attributeValue("type");
                format1 = reader.attributeValue("format");
            }
        }

        // get the value of these fields from the dataset file
        {
            XmlReader reader(dataFilePath);
            if (reader.readNextStartElement())
            {
                root2 = reader.elementName();
                type2 = reader.attributeValue("type");
                format2 = reader.attributeValue("format");
            }
        }

        // return compatible if these fields are nonempty and equal
        // (should be extended some time with more subtle format version matching)
        return !root1.empty() && !type1.empty() && !format1.empty() && root1 == root2 && type1 == type2
               && format1 == format2;
    }

    // catch XML-reader errors and return "incompatible"
    catch (const FatalError&)
    {}
    return false;
}

////////////////////////////////////////////////////////////////////

SchemaDef::SchemaDef(string filePath)
{
    XmlReader reader(filePath);

    // read the root element and its producer attribute
    if (!reader.readNextStartElement() || reader.elementName() != "smile-schema")
        reader.throwError("Expected a 'smile-schema' element");
    _producer = reader.attributeValue("producer");

    // read the Schema element and its attributes
    if (!reader.readNextStartElement() || reader.elementName() != "Schema")
        reader.throwError("Expected a 'Schema' element");
    _name = reader.attributeValue("name");
    _title = reader.attributeValue("title");
    _version = reader.attributeValue("version");
    _extension = reader.attributeValue("extension");
    _root = reader.attributeValue("root");
    _type = reader.attributeValue("type");
    _format = reader.attributeValue("format");
    _url = reader.attributeValue("url");

    // read the highest-level schema property elements ("types", "quantities", ...)
    while (reader.readNextStartElement())
    {
        if (reader.elementName() == "types")
        {
            // read the Type elements
            while (reader.readNextStartElement())
            {
                // get the type name and create a TypeDef entry in our map
                string name = reader.attributeValue("name");
                if (name.empty()) reader.throwError("Type must have a non-empty name");
                TypeDef& def = _allTypes[name];
                def.setName(name);

                // make a special note of concrete types
                if (StringUtils::toBool(reader.attributeValue("concrete")))
                {
                    def.setConcrete();
                    _concreteTypes.push_back(name);
                }

                // get the other basic type attributes
                def.setTitle(reader.attributeValue("title"));
                def.setBase(reader.attributeValue("base"));
                def.setAllowedIf(reader.attributeValue("allowedIf"));
                def.setDisplayedIf(reader.attributeValue("displayedIf"));
                def.setInsert(reader.attributeValue("insert"));
                if (!reader.attributeValue("subPropertyIndex").empty())
                    def.setSubPropertyIndex(StringUtils::toInt(reader.attributeValue("subPropertyIndex")));

                // read the property element, if present
                if (reader.readNextStartElement())
                {
                    if (reader.elementName() != "properties") reader.throwError("Expected a 'properties' element");

                    // read any XxxProperty elements
                    while (reader.readNextStartElement())
                    {
                        // get the basic property attributes and construct a property definition
                        string name = reader.attributeValue("name");
                        if (name.empty()) reader.throwError("Property must have a non-empty name");
                        auto& propDef = def.addPropertyDef(reader.elementName(), name, reader.attributeValue("title"));

                        // add further property info from element attributes
                        propDef.setRelevantIf(reader.attributeValue("relevantIf"));
                        propDef.setDisplayedIf(reader.attributeValue("displayedIf"));
                        propDef.setRequiredIf(reader.attributeValue("requiredIf"));
                        propDef.setInsert(reader.attributeValue("insert"));
                        propDef.setDefaultValue(reader.attributeValue("default"));
                        propDef.setMinValue(reader.attributeValue("min"));
                        propDef.setMaxValue(reader.attributeValue("max"));
                        propDef.setQuantity(reader.attributeValue("quantity"));
                        propDef.setBase(reader.attributeValue("base"));

                        // for enumeration properties, read the elements defining the enumeration values
                        if (propDef.type() == "EnumProperty")
                        {
                            if (reader.readNextStartElement())
                            {
                                if (reader.elementName() != "enumValues")
                                    reader.throwError("Expected an 'enumValues' element");

                                // read any EnumValue elements
                                while (reader.readNextStartElement())
                                {
                                    propDef.addEnumeration(reader.attributeValue("name"),
                                                           reader.attributeValue("title"));
                                    reader.skipCurrentElement();
                                }
                            }
                        }

                        // ignore further contents of the XxxProperty element
                        reader.skipCurrentElement();
                    }

                    // read the Type element end tag
                    reader.skipCurrentElement();
                }
            }
        }
        else if (reader.elementName() == "quantities")
        {
            // read the Quantity elements
            while (reader.readNextStartElement())
            {
                // and remember the corresponding quantity name
                string qty = reader.attributeValue("name");

                // read the units element
                if (reader.readNextStartElement())
                {
                    if (reader.elementName() != "units") reader.throwError("Expected a 'units' element");

                    // read any Unit elements
                    while (reader.readNextStartElement())
                    {
                        string unit = reader.attributeValue("name");
                        double factor = StringUtils::toDouble(reader.attributeValue("factor"));
                        if (!factor) factor = 1.;  // default factor of 1
                        double power = StringUtils::toDouble(reader.attributeValue("power"));
                        if (!power) power = 1.;  // default power of 1
                        double offset = StringUtils::toDouble(reader.attributeValue("offset"));
                        _unitDef.addUnit(qty, unit, factor, power, offset);
                        reader.skipCurrentElement();
                    }
                }

                // read the Quantity element end tag (and ignore any extra information that shouldn't be there anyway)
                reader.skipCurrentElement();
            }
        }
        else if (reader.elementName() == "unitSystems")
        {
            // read the UnitSystem elements
            while (reader.readNextStartElement())
            {
                // and remember the corresponding unit system name
                string unitSystem = reader.attributeValue("name");

                // read the defaultUnits element
                if (reader.readNextStartElement())
                {
                    if (reader.elementName() != "defaultUnits") reader.throwError("Expected a 'defaultUnits' element");

                    // read any DefaultUnit elements
                    while (reader.readNextStartElement())
                    {
                        string quantity = reader.attributeValue("quantity");
                        string unit = reader.attributeValue("unit");
                        _unitDef.addDefaultUnit(unitSystem, quantity, unit);
                        reader.skipCurrentElement();
                    }
                }

                // read the Quantity element end tag (and ignore any extra information that shouldn't be there anyway)
                reader.skipCurrentElement();
            }
        }
        else
            reader.throwError("Expected a 'types', 'quantities', or 'unitSystems' element");
    }

    // verify the end tag of the top-level element
    if (reader.elementName() != "smile-schema") reader.throwError("Expected a 'smile-schema' element end-tag");
}

////////////////////////////////////////////////////////////////////

SchemaDef::SchemaDef(string name, string title, string version, string extension, string root, string type,
                     string format, string url)
    : _name(name), _title(title), _version(version), _extension(extension), _root(root), _type(type), _format(format),
      _url(url)
{}

////////////////////////////////////////////////////////////////////

TypeDef& SchemaDef::addTypeDef(string name, string base, string title, TypeDef::Instantiator instantiator)
{
    // verify that the type is not already there
    if (_allTypes.count(name)) throw FATALERROR("Type '" + name + "' is already defined in the schema");

    // add the type and set basic attributes
    TypeDef& def = _allTypes[name];
    def.setName(name);
    def.setBase(base);
    def.setTitle(title);

    // make a special note of concrete types
    if (instantiator)
    {
        def.setInstantiator(instantiator);
        def.setConcrete();
        _concreteTypes.push_back(name);
    }
    return def;
}

////////////////////////////////////////////////////////////////////

void SchemaDef::loadUnitDef(const UnitDef& unitDef)
{
    _unitDef = unitDef;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // writes the XML stream for the given type definition
    void writeTypeDef(XmlWriter& writer, const TypeDef& tdef)
    {
        // start 'Type' element
        writer.writeStartElement("Type");
        writer.writeAttribute("name", tdef.name());
        writer.writeAttribute("base", tdef.base());
        writer.writeAttribute("title", tdef.title());
        if (tdef.concrete()) writer.writeAttribute("concrete", "true");
        if (tdef.subPropertyIndex() >= 0)
            writer.writeAttribute("subPropertyIndex", std::to_string(tdef.subPropertyIndex()));
        if (!tdef.allowedIf().empty()) writer.writeAttribute("allowedIf", tdef.allowedIf());
        if (!tdef.displayedIf().empty()) writer.writeAttribute("displayedIf", tdef.displayedIf());
        if (!tdef.insert().empty()) writer.writeAttribute("insert", tdef.insert());

        // if this type defines properties
        if (tdef.propertyDefs().size())
        {
            // start 'properties' element,
            writer.writeStartElement("properties");
            writer.writeAttribute("type", "Property");

            // loop over all property definitions
            for (const auto& pdef : tdef.propertyDefs())
            {
                // start Property element
                writer.writeStartElement(pdef.type());
                writer.writeAttribute("name", pdef.name());
                writer.writeAttribute("title", pdef.title());
                if (!pdef.relevantIf().empty()) writer.writeAttribute("relevantIf", pdef.relevantIf());
                if (!pdef.displayedIf().empty()) writer.writeAttribute("displayedIf", pdef.displayedIf());
                if (!pdef.requiredIf().empty()) writer.writeAttribute("requiredIf", pdef.requiredIf());
                if (!pdef.insert().empty()) writer.writeAttribute("insert", pdef.insert());
                if (!pdef.defaultValue().empty()) writer.writeAttribute("default", pdef.defaultValue());
                if (!pdef.minValue().empty()) writer.writeAttribute("min", pdef.minValue());
                if (!pdef.maxValue().empty()) writer.writeAttribute("max", pdef.maxValue());
                if (!pdef.quantity().empty()) writer.writeAttribute("quantity", pdef.quantity());
                if (!pdef.base().empty()) writer.writeAttribute("base", pdef.base());

                // if the property is an enumeration
                if (pdef.enumNames().size())
                {
                    // start 'enumValues' element,
                    writer.writeStartElement("enumValues");
                    writer.writeAttribute("type", "EnumValue");

                    // loop over all property definitions
                    size_t n = pdef.enumNames().size();
                    for (size_t i = 0; i != n; ++i)
                    {
                        writer.writeStartElement("EnumValue");
                        writer.writeAttribute("name", pdef.enumNames()[i]);
                        writer.writeAttribute("title", pdef.enumTitles()[i]);
                        writer.writeEndElement();
                    }

                    // end 'enumValues' element
                    writer.writeEndElement();
                }

                // end Property element
                writer.writeEndElement();
            }

            // end 'properties' element
            writer.writeEndElement();
        }

        // end 'Type' element
        writer.writeEndElement();
    }
}

////////////////////////////////////////////////////////////////////

void SchemaDef::save(string filePath, string producer) const
{
    XmlWriter writer(filePath);

    // write document header
    writer.writeStartDocument();
    writer.writeComment(" A SMILE schema file © Astronomical Observatory, Ghent University ");

    // start root element
    writer.writeStartElement("smile-schema");
    writer.writeAttribute("type", "Schema");
    writer.writeAttribute("format", "1.2");
    writer.writeAttribute("producer", !producer.empty() ? producer : "SMILE");
    writer.writeAttribute("time", System::timestamp(true));

    // start 'Schema' element
    writer.writeStartElement("Schema");
    writer.writeAttribute("name", _name);
    writer.writeAttribute("title", _title);
    writer.writeAttribute("version", _version);
    writer.writeAttribute("extension", _extension);
    writer.writeAttribute("root", _root);
    writer.writeAttribute("type", _type);
    writer.writeAttribute("format", _format);
    writer.writeAttribute("url", _url);

    // ----- types and properties

    // start 'types' element
    writer.writeStartElement("types");
    writer.writeAttribute("type", "Type");

    // loop over all type definitions that are not concrete, in alphabetical order
    writer.writeComment("Non-concrete types, in alphabetical order");
    for (const auto& pair : _allTypes)
    {
        if (!pair.second.concrete()) writeTypeDef(writer, pair.second);
    }

    // loop over all type definitions that are concrete, in order of definition
    writer.writeComment("Concrete types, in order of addition to the registry");
    for (const string& name : _concreteTypes)
    {
        writeTypeDef(writer, _allTypes.at(name));
    }

    // end 'types' element
    writer.writeEndElement();

    // ----- quantities

    // if the schema defines quantities
    if (_unitDef._quantities.size())
    {
        // start 'quantities' element
        writer.writeStartElement("quantities");
        writer.writeAttribute("type", "Quantity");

        // loop over all quantity definitions
        for (const auto& qpair : _unitDef._quantities)
        {
            // start 'Quantity' and 'units' elements
            writer.writeStartElement("Quantity");
            writer.writeAttribute("name", qpair.first);
            writer.writeStartElement("units");
            writer.writeAttribute("type", "Unit");

            // loop over all units
            for (const auto& upair : qpair.second)
            {
                // write "Unit" element
                const string& unitname = upair.first;
                double factor, power, offset;
                std::tie(factor, power, offset) = upair.second;
                writer.writeStartElement("Unit");
                writer.writeAttribute("name", unitname);
                writer.writeAttribute("factor", StringUtils::toString(factor));
                if (power != 1.) writer.writeAttribute("power", StringUtils::toString(power));
                if (offset != 0.) writer.writeAttribute("offset", StringUtils::toString(offset));
                writer.writeEndElement();
            }

            // end 'units' and 'Quantity' elements
            writer.writeEndElement();
            writer.writeEndElement();
        }

        // end 'quantities' element
        writer.writeEndElement();
    }

    // ----- unit systems

    // if the schema defines unit systems
    if (_unitDef._unitSystems.size())
    {
        // start 'unitSystems' element
        writer.writeStartElement("unitSystems");
        writer.writeAttribute("type", "UnitSystem");

        // loop over all unit system definitions
        for (const auto& uspair : _unitDef._unitSystems)
        {
            // start 'UnitSystem' and 'defaultUnits' elements
            writer.writeStartElement("UnitSystem");
            writer.writeAttribute("name", uspair.first);
            writer.writeStartElement("defaultUnits");
            writer.writeAttribute("type", "DefaultUnit");

            // loop over all default units
            for (const auto& dupair : uspair.second)
            {
                // write "DefaultUnit" element
                writer.writeStartElement("DefaultUnit");
                writer.writeAttribute("quantity", dupair.first);
                writer.writeAttribute("unit", dupair.second);
                writer.writeEndElement();
            }

            // end 'units' and 'Quantity' elements
            writer.writeEndElement();
            writer.writeEndElement();
        }

        // end 'unitSystems' element
        writer.writeEndElement();
    }

    // -----

    // end 'Schema' element, root element, and document
    writer.writeEndElement();
    writer.writeEndElement();
    writer.writeEndDocument();
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaProducer() const
{
    return _producer;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaName() const
{
    return _name;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaTitle() const
{
    return _title;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaVersion() const
{
    return _version;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaExtension() const
{
    return _extension;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaRoot() const
{
    return _root;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaType() const
{
    return _type;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaFormat() const
{
    return _format;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::schemaUrl() const
{
    return _url;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::title(string type) const
{
    return typeDef(type).title();
}

////////////////////////////////////////////////////////////////////

vector<string> SchemaDef::titles(const vector<string>& types) const
{
    vector<string> result;
    for (string type : types)
    {
        result.push_back(title(type));
    }
    return result;
}

////////////////////////////////////////////////////////////////////

bool SchemaDef::inherits(string childType, string parentType) const
{
    string type = childType;
    while (!type.empty())
    {
        if (type == parentType) return true;
        type = typeDef(type).base();
    }
    return false;
}

////////////////////////////////////////////////////////////////////

vector<string> SchemaDef::ascendants(string type) const
{
    vector<string> result;
    while (!type.empty())
    {
        result.push_back(type);
        type = typeDef(type).base();
    }
    return result;
}

////////////////////////////////////////////////////////////////////

vector<string> SchemaDef::descendants(string type) const
{
    vector<string> result;
    for (string candidate : _concreteTypes)
    {
        if (inherits(candidate, type)) result.push_back(candidate);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

vector<string> SchemaDef::properties(string type) const
{
    vector<string> result;
    while (!type.empty())
    {
        auto& def = typeDef(type);

        // get the index of the first base property that should be added *after* the sub-properties
        int insertIndex = def.subPropertyIndex();
        if (insertIndex < 0 || insertIndex > def.numSubProperties()) insertIndex = def.numSubProperties();

        // insert the base properties at the front or at the end depending on the above index
        int currentIndex = 0;
        for (auto& propDef : def.propertyDefs())
        {
            if (currentIndex < insertIndex)
                result.insert(result.cbegin() + currentIndex, propDef.name());
            else
                result.push_back(propDef.name());
            currentIndex++;
        }

        type = def.base();
    }
    return result;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::definingType(string type, string property) const
{
    while (!type.empty())
    {
        auto& def = typeDef(type);
        for (auto& propDef : def.propertyDefs())
        {
            if (propDef.name() == property) return type;
        }
        type = def.base();
    }
    throw FATALERROR("Property '" + property + "' is not defined in the schema for this type");
}

////////////////////////////////////////////////////////////////////

string SchemaDef::propertyTitle(string type, string property) const
{
    return propertyDef(type, property).title();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // combine two Boolean expressions into a single one with the "and" operator
    void addAndSegment(string& condition, string segment)
    {
        if (!segment.empty())
        {
            // conditions are evaluated from left to right,
            // so there is no need for parenthesis around the left-hand-side operand
            if (!condition.empty()) condition += "&";
            condition += "(" + segment + ")";
        }
    }
}

////////////////////////////////////////////////////////////////////

string SchemaDef::allowed(string type) const
{
    // combine all "allowedIf" attributes inherited by the type
    string condition;
    while (!type.empty())
    {
        auto& def = typeDef(type);
        addAndSegment(condition, def.allowedIf());
        type = def.base();
    }
    return condition;
}

////////////////////////////////////////////////////////////////////

string SchemaDef::allowedAndDisplayed(string type) const
{
    // combine all "allowedIf" and "displayedIf" attributes inherited by the type
    string condition;
    while (!type.empty())
    {
        auto& def = typeDef(type);
        addAndSegment(condition, def.allowedIf());
        addAndSegment(condition, def.displayedIf());
        type = def.base();
    }
    return condition;
}

////////////////////////////////////////////////////////////////////

vector<string> SchemaDef::toBeInserted(string type) const
{
    vector<string> result;
    while (!type.empty())
    {
        auto& def = typeDef(type);
        if (!def.insert().empty()) result.push_back(def.insert());
        type = def.base();
    }
    return result;
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<Item> SchemaDef::createItem(string type) const
{
    // throws if type is not defined or is abstract
    auto& def = typeDef(type);
    if (!def.concrete()) throw FATALERROR("Can't instantiate abstract type " + type);

    // instantiate real or ghost item
    auto instantiator = def.instantiator();
    if (instantiator)
        return std::unique_ptr<Item>(instantiator());
    else
        return std::make_unique<GhostItem>(type);
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<PropertyHandler> SchemaDef::createPropertyHandler(Item* item, string property,
                                                                  NameManager* nameMgr) const
{
    // get the property definition (throws if not found)
    auto& propDef = propertyDef(item->type(), property);

    // construct handler of subclass corresponding to property type
    string type = propDef.type();
    if (type == "StringProperty") return std::make_unique<StringPropertyHandler>(item, &propDef, this, nameMgr);
    if (type == "BoolProperty") return std::make_unique<BoolPropertyHandler>(item, &propDef, this, nameMgr);
    if (type == "IntProperty") return std::make_unique<IntPropertyHandler>(item, &propDef, this, nameMgr);
    if (type == "EnumProperty") return std::make_unique<EnumPropertyHandler>(item, &propDef, this, nameMgr);
    if (type == "DoubleProperty") return std::make_unique<DoublePropertyHandler>(item, &propDef, this, nameMgr);
    if (type == "DoubleListProperty") return std::make_unique<DoubleListPropertyHandler>(item, &propDef, this, nameMgr);
    if (type == "ItemProperty") return std::make_unique<ItemPropertyHandler>(item, &propDef, this, nameMgr);
    if (type == "ItemListProperty") return std::make_unique<ItemListPropertyHandler>(item, &propDef, this, nameMgr);
    throw FATALERROR("Unsupported property type: " + type);
}

////////////////////////////////////////////////////////////////////

bool SchemaDef::has(string qty) const
{
    return _unitDef.has(qty);
}

////////////////////////////////////////////////////////////////////

bool SchemaDef::has(string qty, string unit) const
{
    return _unitDef.has(qty, unit);
}

////////////////////////////////////////////////////////////////////

double SchemaDef::in(string qty, string unit, double value) const
{
    return _unitDef.in(qty, unit, value);
}

////////////////////////////////////////////////////////////////////

double SchemaDef::out(string qty, string unit, double value) const
{
    return _unitDef.out(qty, unit, value);
}

////////////////////////////////////////////////////////////////////

string SchemaDef::unit(string qty, string unitSystem) const
{
    return _unitDef.unit(qty, unitSystem);
}

////////////////////////////////////////////////////////////////////

bool SchemaDef::hasMultipleUnitSystems() const
{
    return _unitDef._unitSystems.size() > 1;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if both lists contain the same strings in arbitrary order
    // !! assumes that the lists do not contain duplicate strings
    // !! performance does not scale well for long lists
    bool sameStringSets(const vector<string>& va, const vector<string>& vb)
    {
        if (va.size() != vb.size()) return false;
        for (const auto& a : va)
            if (!StringUtils::contains(vb, a)) return false;
        return true;
    }
}

////////////////////////////////////////////////////////////////////

string SchemaDef::unitSystemBase() const
{
    if (_unitDef._unitSystems.size() == 0) throw FATALERROR("No unit system provided in the schema");
    if (_unitDef._unitSystems.size() == 1) return _unitDef._unitSystems.cbegin()->first;

    // build a list of unit system names (there are at least two)
    vector<string> unitSystems;
    for (const auto& pair : _unitDef._unitSystems) unitSystems.push_back(pair.first);

    // consider the ascendants of one of the unit systems as candidates for the common base
    for (const auto& candidate : ascendants(unitSystems[0]))
    {
        if (sameStringSets(descendants(candidate), unitSystems)) return candidate;
    }
    throw FATALERROR("No common base type found for the unit systems in the schema");
}

////////////////////////////////////////////////////////////////////

const TypeDef& SchemaDef::typeDef(string type) const
{
    auto pair = _allTypes.find(type);
    if (pair == _allTypes.cend()) throw FATALERROR("Type '" + type + "' is not defined in the schema");
    return pair->second;
}

////////////////////////////////////////////////////////////////////

const PropertyDef& SchemaDef::propertyDef(string type, string property) const
{
    while (!type.empty())
    {
        auto& def = typeDef(type);
        for (auto& propDef : def.propertyDefs())
        {
            if (propDef.name() == property) return propDef;
        }
        type = def.base();
    }
    throw FATALERROR("Property '" + property + "' is not defined in the schema for this type");
}

////////////////////////////////////////////////////////////////////

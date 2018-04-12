/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XmlHierarchyCreator.hpp"
#include "BoolPropertyHandler.hpp"
#include "DoubleListPropertyHandler.hpp"
#include "DoublePropertyHandler.hpp"
#include "EnumPropertyHandler.hpp"
#include "FatalError.hpp"
#include "IntPropertyHandler.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemPropertyHandler.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SchemaDef.hpp"
#include "StringPropertyHandler.hpp"
#include "StringUtils.hpp"
#include "XmlReader.hpp"
#include <sstream>

////////////////////////////////////////////////////////////////////

namespace
{
    // This function removes enclosing square brackets and the corresponding label from a value string, if present.
    // In other words, it performs the transformation "[label:value]" --> "value".
    // This is necessary to support the mechanism used by FitSKIRT to tag attribute values.
    string removeBrackets(string value)
    {
        if (value.length()>=3 && value.front() == '[' && value.back() == ']' && value.find(':') != string::npos)
        {
            value.erase(value.length()-1);
            value.erase(0, value.find(':')+1);
        }
        return value;
    }

    // ----------------------------------------------------------

    // Forward declarations; see function definitions at the end of this anonymous namespace
    void setupProperties(Item* item, const SchemaDef* schema, XmlReader& reader);
    void setPropertiesToDefaults(Item* item, const SchemaDef* schema, XmlReader& reader);

    // ----------------------------------------------------------

    // The functions in this class are part of the visitor pattern.
    // They set the value of a property read from the current position in an XML stream.
    class ReaderPropertySetter : public PropertyHandlerVisitor
    {
    private:
        XmlReader& _reader;

    public:
        ReaderPropertySetter(XmlReader& reader)
            : _reader(reader) { }

        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            string value = _reader.attributeValue(handler->name());
            if (!value.empty()) handler->setValue(value);
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override
        {
            string value = removeBrackets(_reader.attributeValue(handler->name()));
            if (StringUtils::isValidBool(value))
            {
                handler->setValue(StringUtils::toBool(value));
            }
            else _reader.throwError("Value '" + value + "' for property '" + handler->name()
                                    + "' can't be converted to bool");
        }

        void visitPropertyHandler(IntPropertyHandler* handler) override
        {
            string value = removeBrackets(_reader.attributeValue(handler->name()));
            if (StringUtils::isValidInt(value))
            {
                int ivalue = StringUtils::toInt(value);
                if (ivalue >= handler->minValue() && ivalue <= handler->maxValue())
                {
                    handler->setValue(ivalue);
                }
                else _reader.throwError("Value '" + value + "' for property '" + handler->name()
                                       + "' is out of range");
            }
            else _reader.throwError("Value '" + value + "' for property '" + handler->name()
                                    + "' can't be converted to integer");
        }

        void visitPropertyHandler(EnumPropertyHandler* handler) override
        {
            string value = removeBrackets(_reader.attributeValue(handler->name()));
            if (handler->isValid(value))
            {
                handler->setValue(value);
            }
            else _reader.throwError("Value '" + value + "' for property '" + handler->name()
                                    + "' is an invalid enumeration key");
        }

        void visitPropertyHandler(DoublePropertyHandler* handler) override
        {
            string value = removeBrackets(_reader.attributeValue(handler->name()));
            if (handler->isValidDouble(value))
            {
                double dvalue = handler->toDouble(value);
                if (handler->isInRange(dvalue))
                {
                    handler->setValue(dvalue);
                }
                else _reader.throwError("Value '" + value + "' for property '" + handler->name()
                                       + "' is out of range " + handler->rangeDescription());
            }
            else _reader.throwError("Value '" + value + "' for property '" + handler->name()
                                    + "' can't be converted to double");
        }

        void visitPropertyHandler(DoubleListPropertyHandler* handler) override
        {
            string value = removeBrackets(_reader.attributeValue(handler->name()));
            if (handler->isValidDoubleList(value))
            {
                auto lvalue = handler->toDoubleList(value);
                if (handler->isInRange(lvalue))
                {
                    handler->setValue(lvalue);
                }
                else _reader.throwError("Value(s) '" + value + "' for property '" + handler->name()
                                       + "' is (are) out of range " + handler->rangeDescription());
            }
            else _reader.throwError("Value '" + value + "' for property '" + handler->name()
                                    + "' can't be converted to list of doubles");
        }

        void visitPropertyHandler(ItemPropertyHandler* handler) override
        {
            // verify the property base type
            string baseType = _reader.attributeValue("type");
            if (baseType != handler->baseType())
                _reader.throwError("Type '" + baseType + "' does not match base type '" + handler->baseType()
                                        + "' for property '" + handler->name() + "'");

            // read the element for the item and verify its actual type
            if (!_reader.readNextStartElement())
                _reader.throwError("Expected element for item with type inheriting '" + baseType + "'");
            string type = _reader.elementName();
            if (!handler->schema()->inherits(type, baseType))
                _reader.throwError("Item with type '" + type + "' does not inherit '" + baseType + "'");

            // create the item and set the property
            bool success = handler->setToNewItemOfType(type);
            if (!success) _reader.throwError("Can't create item of type " + type);

            // recursively handle the newly created item
            setupProperties(handler->value(), handler->schema(), _reader);

            // process the end of the property element
            if (_reader.readNextStartElement())
                _reader.throwError("Unexpected element '" + _reader.elementName() + "'");
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            // verify the property base type
            string baseType = _reader.attributeValue("type");
            if (baseType != handler->baseType())
                _reader.throwError("Type '" + baseType + "' does not match base type '" + handler->baseType()
                                        + "' for property '" + handler->name() + "'");

            // initialize the property value
            handler->setToEmpty();

            // process all elements at this level, one by one
            while (_reader.readNextStartElement())
            {
                // verify the actual type of the item
                string type = _reader.elementName();
                if (!handler->schema()->inherits(type, baseType))
                    _reader.throwError("Item with type '" + type + "' does not inherit '" + baseType + "'");

                // create the item and add it to the property
                bool success = handler->addNewItemOfType(type);
                if (!success) _reader.throwError("Can't create item of type " + type);

                // recursively handle the newly created item
                setupProperties(handler->value().back(), handler->schema(), _reader);
            }
        }
    };

    // ----------------------------------------------------------

    // The functions in this class are part of the visitor pattern.
    // They set the value of a property to its default value (assuming it has one).
    class DefaultPropertySetter : public PropertyHandlerVisitor
    {
    private:
        XmlReader& _reader;

    public:
        DefaultPropertySetter(XmlReader& reader)
            : _reader(reader) { }

        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(IntPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(EnumPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(DoublePropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(DoubleListPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(ItemPropertyHandler* handler) override
        {
            if (!handler->isOptional() && handler->isRelevant())
            {
                // create the default item and set the property
                bool success = handler->setToNewItemOfType(handler->defaultType());
                if (!success) _reader.throwError("Can't create item of type " + handler->defaultType());

                // recursively default-construct the properties of the new item
                setPropertiesToDefaults(handler->value(), handler->schema(), _reader);
            }
            else handler->setToNull();
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            handler->setToEmpty();
            if (!handler->isOptional() && handler->isRelevant())
            {
                // create the default item and add it to the property
                bool success = handler->addNewItemOfType(handler->defaultType());
                if (!success) _reader.throwError("Can't create item of type " + handler->defaultType());

                // recursively default-construct the properties of the new item
                setPropertiesToDefaults(handler->value().back(), handler->schema(), _reader);
            }
        }
    };

    // ----------------------------------------------------------

    // This function recursively sets up the properties of the specified SMILE data item and its children.
    // Setting up the properties at the current level poses a dual challenge:
    //  - we need to process scalar and compound properties in the order provided in the XML file
    //  - we need to ensure that all properties defined for the item are handled
    // To achieve this, we proceed as follows:
    //  - build a dictionary of handlers for all defined properties
    //  - process the properties in the XML file
    //  - the handlers keep track of whether a value has been set
    //  - process the "virgin" handlers to set default values or complain if there is no default
    // Actually setting the values is accomplished by asking each of the handlers to accept an appropriate
    // PropertyHandlerVisitor instance as a visitor, which causes a call-back to the visitPropertyHandler()
    // function with the corresponding PropertyHandler type.
    void setupProperties(Item* item, const SchemaDef* schema, XmlReader& reader)
    {
        ReaderPropertySetter readerSetter(reader);
        DefaultPropertySetter defaultSetter(reader);

        // build a dictionary of handlers for all defined properties
        std::unordered_map<string,std::unique_ptr<PropertyHandler>> handlers;
        for (const string& name : schema->properties(item->type()))
        {
            handlers[name] = schema->createPropertyHandler(item, name);
        }

        // process scalar properties (derived from XML attributes)
        for (const string& name : reader.attributeNames())
        {
            if (!handlers.count(name))
                reader.throwError("Item of type '" + item->type() + "' has no property named '" + name + "'");
            auto& handler = handlers[name];
            if (handler->isCompound())
                reader.throwError("Property '" + handler->name() +
                                  "' has a non-compound data type and is given as an xml element");
            handler->acceptVisitor(&readerSetter);
        }

        // process compound properties (derived from XML child elements)
        while (reader.readNextStartElement())
        {
            string name = reader.elementName();
            if (!handlers.count(name))
                reader.throwError("Item of type '" + item->type() + "' has no property named '" + name + "'");
            auto& handler = handlers[name];
            if (!handler->isCompound())
                reader.throwError("Property '" + handler->name() +
                                   "' has a compund data type and is given as an xml attribute");
            handler->acceptVisitor(&readerSetter);
        }

        // honor default property values, in order of schema definition so that relevancy can be determined
        for (const string& name : schema->properties(item->type()))
        {
            auto& handler = handlers[name];
            if (!handler->hasChanged())
            {
                if (!handler->isOptional() && !handler->hasDefaultValue() && handler->isRelevant())
                    reader.throwError("Value for required property '" + handler->name() + "' in item of type '"
                                      + item->type() + "' is not specified and has no default value");
                handler->acceptVisitor(&defaultSetter);
            }
        }
    }

    // ----------------------------------------------------------

    // This function recursively sets all properties of the specified SMILE data item and its children
    // to their default values, without reading anything from the XML stream. The function throws an error
    // if not all properties (recursively) have a default value.
    void setPropertiesToDefaults(Item* item, const SchemaDef* schema, XmlReader& reader)
    {
        DefaultPropertySetter defaultSetter(reader);

        // process properties in order of schema definition so that relevancy can be determined
        for (const string& name : schema->properties(item->type()))
        {
            auto handler = schema->createPropertyHandler(item, name);
            if (!handler->isOptional() && !handler->hasDefaultValue() && handler->isRelevant())
                reader.throwError("Value for required property '" + handler->name() + "' in item of type '"
                                  + item->type() + "' is not specified and has no default value");
            handler->acceptVisitor(&defaultSetter);
        }
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    std::unique_ptr<Item> read(const SchemaDef* schema, XmlReader& reader)
    {
        // read the root element and verify the top-level base type
        if (!reader.readNextStartElement())
            reader.throwError("Can't find XML root element");
        if (reader.elementName() != schema->schemaRoot())
            reader.throwError("Root element is '" + reader.elementName() + "' rather than '" + schema->schemaRoot() + "'");
        string baseType = schema->schemaType();
        if (reader.attributeValue("type") != baseType)
            reader.throwError("Top-level base type '" + reader.attributeValue("type") +
                              "' does not match '" + baseType + "'");

        // read the element for the top-level item and verify its actual type
        if (!reader.readNextStartElement())
            reader.throwError("Expected element for top-level item with type inheriting '" + baseType + "'");
        string type = reader.elementName();
        if (!schema->inherits(type, baseType))
            reader.throwError("Top-level item with type '" + type + "' does not inherit '" + baseType + "'");

        // create the top-level item
        auto topLevelItem = schema->createItem(type);

        // recursively setup all properties of the top-level item and its children
        setupProperties(topLevelItem.get(), schema, reader);

        // process the end of the root element
        if (reader.readNextStartElement())
            reader.throwError("Unexpected element '" + reader.elementName() + "'");

        return topLevelItem;
    }
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<Item> XmlHierarchyCreator::readFile(const SchemaDef* schema, string filepath)
{
    // construct the XML reader and call the common read() function
    XmlReader reader(filepath);
    return read(schema, reader);
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<Item> XmlHierarchyCreator::readString(const SchemaDef* schema, string contents, string description)
{
    // construct the XML reader and call the common read() function
    std::istringstream stream(contents);
    XmlReader reader(stream, description);
    return read(schema, reader);
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XmlHierarchyWriter.hpp"
#include "BoolPropertyHandler.hpp"
#include "DoubleListPropertyHandler.hpp"
#include "DoublePropertyHandler.hpp"
#include "EnumPropertyHandler.hpp"
#include "FatalError.hpp"
#include "IntPropertyHandler.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemPropertyHandler.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SchemaDef.hpp"
#include "StringPropertyHandler.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "XmlWriter.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // Forward declaration; see function definition at the end of this anonymous namespace
    void writeProperties(Item* item, const SchemaDef* schema, XmlWriter& writer);

    // ----------------------------------------------------------

    // The functions in this class are part of the visitor pattern initiated by the setupProperties() function.
    // They write the appropriate XML information for the specified property.
    class PropertyWriter : public PropertyHandlerVisitor
    {
    private:
        const SchemaDef* _schema;
        XmlWriter& _writer;

    public:
        PropertyWriter(const SchemaDef* schema, XmlWriter& writer) : _schema(schema), _writer(writer) {}

        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            _writer.writeAttribute(handler->name(), handler->value());
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override
        {
            _writer.writeAttribute(handler->name(), StringUtils::toString(handler->value()));
        }

        void visitPropertyHandler(IntPropertyHandler* handler) override
        {
            _writer.writeAttribute(handler->name(), StringUtils::toString(handler->value()));
        }

        void visitPropertyHandler(EnumPropertyHandler* handler) override
        {
            _writer.writeAttribute(handler->name(), handler->value());
        }

        void visitPropertyHandler(DoublePropertyHandler* handler) override
        {
            _writer.writeAttribute(handler->name(), handler->toString(handler->value()));
        }

        void visitPropertyHandler(DoubleListPropertyHandler* handler) override
        {
            _writer.writeAttribute(handler->name(), handler->toString(handler->value()));
        }

        void visitPropertyHandler(ItemPropertyHandler* handler) override
        {
            if (handler->value())
            {
                // start an element for the property
                _writer.writeStartElement(handler->name());
                _writer.writeAttribute("type", handler->baseType());

                // handle the item pointed to by the property
                writeProperties(handler->value(), _schema, _writer);

                // end the element for the property
                _writer.writeEndElement();
            }
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            auto items = handler->value();
            if (!items.empty())
            {
                // start an element for the property
                _writer.writeStartElement(handler->name());
                _writer.writeAttribute("type", handler->baseType());

                // handle the items pointed to by the property
                for (Item* item : items)
                {
                    writeProperties(item, _schema, _writer);
                }

                // end the element for the property
                _writer.writeEndElement();
            }
        }
    };

    // ----------------------------------------------------------

    // This function recursively writes the properties of the specified item and its children.
    // This is accomplished by asking each of the properties to accept an appropriate PropertyHandlerVisitor
    // instance as a visitor, which causes a call-back to the visitPropertyHandler() function with
    // the corresponding PropertyHandler type.
    void writeProperties(Item* item, const SchemaDef* schema, XmlWriter& writer)
    {
        PropertyWriter propertyWriter(schema, writer);

        // start an element for the item
        writer.writeStartElement(item->type());

        // handle the item's properties by distributing them to write methods depending on property type (visitor pattern)
        //  - first handle properties that are written to XML attributes
        //  - then handle properties that are written to XML as child elements
        for (const string& property : schema->properties(item->type()))
        {
            auto handler = schema->createPropertyHandler(item, property, nullptr);
            if (!handler->isCompound()) handler->acceptVisitor(&propertyWriter);
        }
        for (const string& property : schema->properties(item->type()))
        {
            auto handler = schema->createPropertyHandler(item, property, nullptr);
            if (handler->isCompound()) handler->acceptVisitor(&propertyWriter);
        }

        // end the element for the item
        writer.writeEndElement();
    }
}

////////////////////////////////////////////////////////////////////

void XmlHierarchyWriter::write(Item* item, const SchemaDef* schema, string filePath, string producer)
{
    // setup the XML writer and cache some pointers for use in other member functions
    XmlWriter writer(filePath);

    // verify the type of the top-level item
    if (!schema->inherits(item->type(), schema->schemaType()))
        throw FATALERROR("Top-level item type " + item->type() + " does not derive from schema type "
                         + schema->schemaType());

    // write document header and start root element
    writer.writeStartDocument();
    writer.writeComment(" " + StringUtils::toUpperFirst(schema->schemaTitle())
                        + " © Astronomical Observatory, Ghent University ");
    writer.writeStartElement(schema->schemaRoot());
    writer.writeAttribute("type", schema->schemaType());
    writer.writeAttribute("format", schema->schemaFormat());
    writer.writeAttribute("producer", !producer.empty() ? producer : "SMILE");
    writer.writeAttribute("time", System::timestamp(true));

    // recursively write all properties of the top-level item and its children
    writeProperties(item, schema, writer);

    // end root element and document
    writer.writeEndElement();
    writer.writeEndDocument();
}

////////////////////////////////////////////////////////////////////

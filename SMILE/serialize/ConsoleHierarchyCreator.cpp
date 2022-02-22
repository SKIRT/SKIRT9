/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ConsoleHierarchyCreator.hpp"
#include "BoolPropertyHandler.hpp"
#include "Console.hpp"
#include "DoubleListPropertyHandler.hpp"
#include "DoublePropertyHandler.hpp"
#include "EnumPropertyHandler.hpp"
#include "FatalError.hpp"
#include "IntPropertyHandler.hpp"
#include "Item.hpp"
#include "ItemListPropertyHandler.hpp"
#include "ItemPropertyHandler.hpp"
#include "NameManager.hpp"
#include "PropertyDef.hpp"
#include "PropertyHandlerVisitor.hpp"
#include "SchemaDef.hpp"
#include "StringPropertyHandler.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // Forward declarations; see function definitions at the end of this anonymous namespace
    void setPropertiesToDefaults(Item* item, const SchemaDef* schema, NameManager* nameMgr);
    void setupProperties(Item* item, const SchemaDef* schema, NameManager* nameMgr);

    // ----------------------------------------------------------

    // The functions in this class are part of the visitor pattern initiated by the setupProperties() function.
    // They set the value of a property by asking the user at the console for a value.
    class ConsolePropertySetter : public PropertyHandlerVisitor
    {
    public:
        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            string value =
                Console::promptForString("Enter " + handler->title(),
                                         handler->hasDefaultValue() || !handler->isRequired(), handler->defaultValue());
            handler->setValue(value);
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override
        {
            bool value = Console::promptForBool("Do you want to " + handler->title() + "?", handler->hasDefaultValue(),
                                                handler->defaultValue());
            handler->setValue(value);
        }

        void visitPropertyHandler(IntPropertyHandler* handler) override
        {
            int value = Console::promptForInt("Enter " + handler->title(), handler->minValue(), handler->maxValue(),
                                              handler->hasDefaultValue(), handler->defaultValue());
            handler->setValue(value);
        }

        void visitPropertyHandler(EnumPropertyHandler* handler) override
        {
            vector<string> names = handler->values();
            int choice =
                Console::promptForChoice(handler->title(), handler->titlesForValues(), handler->hasDefaultValue(),
                                         StringUtils::indexOf(names, handler->defaultValue()));
            handler->setValue(names[choice]);
        }

        void visitPropertyHandler(DoublePropertyHandler* handler) override
        {
            handler->setValue(ConsoleHierarchyCreator::promptForDouble("Enter ", handler));
        }

        void visitPropertyHandler(DoubleListPropertyHandler* handler) override
        {
            handler->setValue(ConsoleHierarchyCreator::promptForDoubleList("Enter ", handler));
        }

        void visitPropertyHandler(ItemPropertyHandler* handler) override
        {
            // make the user select the appropriate subclass for this property
            vector<string> choices = handler->allowedAndDisplayedDescendants();
            int choice = Console::promptForChoice(handler->title(), handler->schema()->titles(choices),
                                                  handler->hasDefaultValue(),
                                                  StringUtils::indexOf(choices, handler->defaultType()),
                                                  !handler->isRequired(), "or zero to select none");

            if (choice < 0)
            {
                // set the property value to null (i.e. absent) if so requested
                handler->setToNull();
            }
            else
            {
                // create the item and set the property
                bool success = handler->setToNewItemOfType(choices[choice]);
                if (!success) throw FATALERROR("Can't create simulation item of type " + choices[choice]);

                // recursively handle the newly created simulation item
                setupProperties(handler->value(), handler->schema(), handler->nameManager());
            }
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            // get a list of the subclasses for this property
            vector<string> choices = handler->allowedAndDisplayedDescendants();

            // initialize the property value
            handler->setToEmpty();

            // loop to create list of new items
            for (int count = 1;; count++)
            {
                // make the user select the appropriate subclass for this item
                int choice =
                    Console::promptForChoice("item #" + std::to_string(count) + " in " + handler->title() + " list",
                                             handler->schema()->titles(choices), handler->hasDefaultValue(),
                                             StringUtils::indexOf(choices, handler->defaultType()),
                                             count != 1 || !handler->isRequired(), "or zero to terminate the list");

                // terminate the list if requested
                if (choice < 0) return;

                // create the item and set the property
                bool success = handler->addNewItemOfType(choices[choice]);
                if (!success) throw FATALERROR("Can't create simulation item of type " + choices[choice]);

                // recursively handle the newly created simulation item
                setupProperties(handler->value().back(), handler->schema(), handler->nameManager());
            }
        }
    };

    // ----------------------------------------------------------

    // The functions in this class are part of the visitor pattern.
    // They set the value of a property to its default value (assuming it has one).
    class DefaultPropertySetter : public PropertyHandlerVisitor
    {
    public:
        void visitPropertyHandler(StringPropertyHandler* handler) override
        {
            handler->setValue(handler->defaultValue());
        }

        void visitPropertyHandler(BoolPropertyHandler* handler) override { handler->setValue(handler->defaultValue()); }

        void visitPropertyHandler(IntPropertyHandler* handler) override { handler->setValue(handler->defaultValue()); }

        void visitPropertyHandler(EnumPropertyHandler* handler) override { handler->setValue(handler->defaultValue()); }

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
            if (handler->isRelevant() && handler->isRequired())
            {
                // create the default item and set the property
                bool success = handler->setToNewItemOfType(handler->defaultType());
                if (!success) throw FATALERROR("Can't create item of type " + handler->defaultType());

                // recursively default-construct the properties of the new item
                setPropertiesToDefaults(handler->value(), handler->schema(), handler->nameManager());
            }
            else
                handler->setToNull();
        }

        void visitPropertyHandler(ItemListPropertyHandler* handler) override
        {
            handler->setToEmpty();
            if (handler->isRelevant() && handler->isRequired())
            {
                // create the default item and add it to the property
                bool success = handler->addNewItemOfType(handler->defaultType());
                if (!success) throw FATALERROR("Can't create item of type " + handler->defaultType());

                // recursively default-construct the properties of the new item
                setPropertiesToDefaults(handler->value().back(), handler->schema(), handler->nameManager());
            }
        }
    };

    // ----------------------------------------------------------

    // This function recursively sets all properties of the specified SMILE data item and its children
    // to their default values, without reading anything from the console. The function throws an error
    // if not all properties (recursively) have a default value.
    void setPropertiesToDefaults(Item* item, const SchemaDef* schema, NameManager* nameMgr)
    {
        DefaultPropertySetter defaultSetter;

        // process properties in order of schema definition so that relevancy can be determined
        for (const string& name : schema->properties(item->type()))
        {
            auto handler = schema->createPropertyHandler(item, name, nameMgr);
            if (handler->isRelevant() && handler->isRequired() && !handler->hasDefaultValue())
                throw FATALERROR("Value for required property '" + handler->name() + "' in item of type '"
                                 + item->type() + "' is not specified and has no default value");
            handler->acceptVisitor(&defaultSetter);
        }
    }

    // ----------------------------------------------------------

    // This function recursively sets up the properties of the specified SMILE data item and its children.
    // This is accomplished by asking each of the properties to accept an appropriate PropertyHandlerVisitor
    // instance as a visitor, which causes a call-back to the visitPropertyHandler() function with
    // the corresponding PropertyHandler type.
    void setupProperties(Item* item, const SchemaDef* schema, NameManager* nameMgr)
    {
        ConsolePropertySetter consoleSetter;
        DefaultPropertySetter defaultSetter;

        // setup all properties of the item, using a fresh local name space
        nameMgr->pushLocal();
        for (const string& property : schema->properties(item->type()))
        {
            auto handler = schema->createPropertyHandler(item, property, nameMgr);

            // distribute to setup methods depending on property type (using visitor pattern)
            if (handler->isSilent())
                handler->acceptVisitor(&defaultSetter);
            else
                handler->acceptVisitor(&consoleSetter);
        }
        nameMgr->popLocal();
    }
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<Item> ConsoleHierarchyCreator::create(const SchemaDef* schema)
{
    // make the user select a subtype of the root base type
    string rootBaseType = schema->schemaType();
    vector<string> choices = schema->descendants(rootBaseType);
    int choice = Console::promptForChoice(schema->title(rootBaseType), schema->titles(choices));
    string rootType = choices[choice];

    // create the root item
    auto rootItem = schema->createItem(rootType);

    // construct the name manager for this session, and insert the appropriate names for the root item type
    NameManager nameMgr;
    nameMgr.insert(schema->ascendants(rootType));
    nameMgr.insertFromConditionalValue(schema->toBeInserted(rootType));

    // recursively setup all properties of the root item and its children
    setupProperties(rootItem.get(), schema, &nameMgr);
    return rootItem;
}

////////////////////////////////////////////////////////////////////

double ConsoleHierarchyCreator::promptForDouble(string prefix, const DoublePropertyHandler* handler)
{
    // get default and min/max values and verify that default is in range
    bool hasDef = handler->hasDefaultValue();
    double def = handler->defaultValue();
    if (hasDef && !handler->isInRange(def)) throw FATALERROR("Default double value out of range");

    // add prefix, title, hints on min/max and default values to the message
    string message = prefix + handler->title() + " " + handler->rangeDescription();
    if (hasDef) message += " (" + handler->toString(def) + ")";

    while (true)
    {
        // get the input string and attempt conversion
        string input = System::prompt(message);
        double result = handler->toDouble(input);

        // if provided, use the default value instead of an empty string
        if (input.empty() && hasDef) return def;

        // reject conversion errors or empty strings without default
        if (!handler->isValidDouble(input))
        {
            Console::error("Enter a valid floating point number, optionally followed by a space and a unit string");
        }

        // reject out-of-range values
        else if (!handler->isInRange(result))
        {
            Console::error("Enter a number within range " + handler->rangeDescription());
        }

        // return successful result
        else
        {
            return result;
        }
    }
}

////////////////////////////////////////////////////////////////////

vector<double> ConsoleHierarchyCreator::promptForDoubleList(string prefix, const DoubleListPropertyHandler* handler)
{
    // get default and min/max values and verify that default is in range
    bool hasDef = handler->hasDefaultValue() || !handler->isRequired();
    vector<double> def = handler->defaultValue();
    if (hasDef && !handler->isInRange(def)) throw FATALERROR("Default double list value out of range");

    // add prefix, title, hints on min/max and default values to the message
    string message = prefix + handler->title() + " " + handler->rangeDescription();
    if (hasDef) message += " (" + handler->toString(def) + ")";

    while (true)
    {
        // get the input string
        string input = System::prompt(message);

        // if provided, use the default value instead of an empty string
        if (input.empty() && hasDef) return def;

        // reject conversion errors or empty strings without default
        if (!handler->isValidDoubleList(input))
        {
            Console::error("Enter a comma-separated list of floating point numbers, "
                           "each optionally followed by a space and a unit string");
        }

        // reject out-of-range values
        else
        {
            // convert the list of numbers
            vector<double> result = handler->toDoubleList(input);

            // verify that all numbers are within range
            if (!handler->isInRange(result))
            {
                Console::error("Enter numbers within range " + handler->rangeDescription());
            }

            // return successful result
            else
            {
                return result;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

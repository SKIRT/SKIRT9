/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOOLPROPERTYHANDLER_HPP
#define BOOLPROPERTYHANDLER_HPP

#include "PropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties of Boolean types. */
class BoolPropertyHandler : public PropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using PropertyHandler::PropertyHandler;

    // ================== Overriding base class functions ==================

public:
    /** Returns true if the given string can be successfully converted to a value of the property's
        type. For Boolean properties, the function returns true if the string is recognized as
        representing a Boolean by the StringUtils::isValidBool() function, and false otherwise. */
    bool isValidValue(string value) const override;

    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. For Boolean
        properties, the function inserts the target property's name if the current property value
        is true, and does not insert any names if the value is false. In addition, the function
        inserts the names provided in the conditional expression of the "insert" attribute of the
        target property, if any. */
    void insertNames() override;

    /** Accepts the specified property handler visitor. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Specific functions for this property type ==================

public:
    /** Returns the default value for the handled property, or false if unavailable. */
    bool defaultValue() const;

    /** Returns the value of the handled property in the target item. */
    bool value() const;

    /** Sets the value of the handled property in the target item. */
    void setValue(bool value);
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STRINGPROPERTYHANDLER_HPP
#define STRINGPROPERTYHANDLER_HPP

#include "PropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties of string types. */
class StringPropertyHandler : public PropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using PropertyHandler::PropertyHandler;

    // ================== Overriding base class functions ==================

public:
    /** Returns true if the value of the handled property is a non-empty string, or false if the
        value is the empty string. */
    bool isTrueInCondition() const override;

    /** Returns true if the handled property is optional (i.e. it may have an empty value), or
        false if it is required. */
    bool isOptional() const override;

    /** Returns true if the handled property has a valid default value, or false if not. */
    bool hasDefaultValue() const override;

    /** Accepts the specified property handler visitor. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Functionality for this property type ==================

public:
    /** Returns the default value for the handled property, or the empty string if unavailable. */
    string defaultValue() const;

    /** Returns the value of the handled property in the target item. */
    string value() const;

    /** Sets the value of the handled property in the target item. */
    void setValue(string value);
};

////////////////////////////////////////////////////////////////////

#endif

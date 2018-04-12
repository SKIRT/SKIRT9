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
    /** Returns the value of the handled property. */
    bool isTrueInCondition() const override;

    /** Returns true if the handled property has a valid default value, or false if not. */
    bool hasDefaultValue() const override;

    /** Accepts the specified property handler visitor. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Functionality for this property type ==================

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

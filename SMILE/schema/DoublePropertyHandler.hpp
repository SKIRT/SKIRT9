/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DOUBLEPROPERTYHANDLER_HPP
#define DOUBLEPROPERTYHANDLER_HPP

#include "AbstractDoublePropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties that hold a single floating point value with an
    optional unit specification. */
class DoublePropertyHandler : public AbstractDoublePropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using AbstractDoublePropertyHandler::AbstractDoublePropertyHandler;

    // ================== Overriding base class functions ==================

public:
    /** Returns true if the given string can be successfully converted to a value of the property's
        type. For double properties, the function returns true if the string conforms to the syntax
        recognized by the AbstractDoublePropertyHandler::isValidDouble() function, and false
        otherwise. */
    bool isValidValue(string value) const override;

    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. For double
        properties, the function inserts the target property's name if the current property value
        is nonzero, and does not insert any names if the value is zero. In addition, the function
        inserts the names provided in the conditional expression of the "insert" attribute of the
        target property, if any. */
    void insertNames() override;

    /** Accepts the specified visitor. This function is part of the "visitor" design pattern
        implementation used to handle properties of various types. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Specific functions for this property type ==================

public:
    /** Returns the default value for the handled property, or zero if unavailable. */
    double defaultValue() const;

    /** Returns the value of the handled property in the target item. */
    double value() const;

    /** Sets the value of the handled property in the target item. */
    void setValue(double value);
};

////////////////////////////////////////////////////////////////////

#endif

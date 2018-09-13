/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INTPROPERTYHANDLER_HPP
#define INTPROPERTYHANDLER_HPP

#include "PropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties of integer types. */
class IntPropertyHandler : public PropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using PropertyHandler::PropertyHandler;

    // ================== Overriding base class functions ==================

public:
    /** Returns true if the given string can be successfully converted to a value of the property's
        type. For integer properties, the function returns true if the string conforms to the
        regular syntax for a decimal integer, and false otherwise. */
    bool isValidValue(string value) const override;

    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. For integer
        properties, the function inserts the target property's name if the current property value
        is nonzero, and does not insert any names if the value is zero. In addition, the function
        inserts the names provided in the conditional expression of the "insert" attribute of the
        target property, if any. */
    void insertNames() override;

    /** Accepts the specified property handler visitor. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Specific functions for this property type ==================

public:
    /** Returns the default value for the handled property, or zero if unavailable. */
    int defaultValue() const;

    /** Returns the minimum value for the handled property. If no minimum value is specified in the
        property definition, the function returns a default minimum value close to the smallest
        representable integer (i.e. a negative integer with a large absolute value). */
    int minValue() const;

    /** Returns the maximum value for the handled property. If no maximum value is specified in the
        property definition, the function returns a default maximum value close to the largest
        representable integer. */
    int maxValue() const;

    /** Returns the value of the handled property in the target item. */
    int value() const;

    /** Sets the value of the handled property in the target item. */
    void setValue(int value);
};

////////////////////////////////////////////////////////////////////

#endif

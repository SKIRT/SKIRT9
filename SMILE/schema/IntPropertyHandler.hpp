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
    /** Returns true if the value of the handled property is nonzero, and false if it is zero. */
    bool isTrueInCondition() const override;

    /** Returns true if the handled property has a valid default value, or false if not. */
    bool hasDefaultValue() const override;

    /** Accepts the specified property handler visitor. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Functionality for this property type ==================

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

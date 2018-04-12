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
    /** Returns true if the handled property has a valid default value, or false if not. */
    bool hasDefaultValue() const override;

    /** Accepts the specified visitor. This function is part of the "visitor" design pattern
        implementation used to handle properties of various types. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Functionality for this property type ==================

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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DOUBLELISTPROPERTYHANDLER_HPP
#define DOUBLELISTPROPERTYHANDLER_HPP

#include "AbstractDoublePropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties that hold a list of floating point values with an
    optional unit specification, seperated by commas. */
class DoubleListPropertyHandler : public AbstractDoublePropertyHandler
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
    vector<double> defaultValue() const;

    /** Returns the value of the handled property in the target item. */
    vector<double> value() const;

    /** Sets the value of the handled property in the target item. */
    void setValue(vector<double> value);
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSTRACTITEMPROPERTYHANDLER_HPP
#define ABSTRACTITEMPROPERTYHANDLER_HPP

#include "PropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** AbstractItemPropertyHandler serves as an abstract base class for handling SMILE data item
    properties that hold one or more pointers to other SMILE data items. */
class AbstractItemPropertyHandler : public PropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using PropertyHandler::PropertyHandler;

    // ================== Overriding base class functions ==================

public:
    /** Returns true if the given string can be successfully converted to a value of the property's
        type. For item and item list properties, the function returns true if the string matches
        the type name of one of the currently allowed descendants, and false otherwise. */
    bool isValidValue(string value) const override;

    /** Returns true, indicating that the handled property type is compound, i.e. it aggregates
        other items that are part of the item hierarchy. */
    bool isCompound() const override;

    // ================== Specific functions for this property type ==================

public:
    /** Returns the base type of the item being pointed to by the handled property. */
    string baseType() const;

    /** Returns the default type for the handled property, or empty if unavailable. */
    string defaultType() const;

    /** Returns a list of item types, in order of addition to the item registry, which inherit the
        base type of the target item and which are allowed and should be displayed for the current
        dataset configuration. For each candidate type, the function obtains a compound conditional
        expression by concatenating the "allowedIf" and "displayedIf" attribute values for that
        type and for any of its base types, recursively, with the AND operator into a single
        Boolean expression. After replacing names currently present in the global or local name
        sets by true, and other names by false, it then evaluates this expression to determine
        whether the type is to be included in the returned list. */
    vector<string> allowedAndDisplayedDescendants();
};

////////////////////////////////////////////////////////////////////

#endif

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
    /** Returns true if the handled property is optional (i.e. its value may be a null pointer or
        an empty list), or false if not. */
    bool isOptional() const override;

    /** Returns true if the handled property has a valid default value, or false if not. */
    bool hasDefaultValue() const override;

    /** Returns true, indicating that the handled property type is compound, i.e. it aggregates
        other items that are part of the item hierarchy. */
    bool isCompound() const override;

    // ================== Functionality for this property type ==================

public:
    /** Returns the base type of the item being pointed to by the handled property. */
    string baseType() const;

    /** Returns the default type for the handled property, or empty if unavailable. */
    string defaultType() const;

    // ================== Utilities for editing ==================

public:
    /** This function returns a list of item types, in order of addition to the item registry,
        which inherit the base type of the target item and which are allowed according to
        conditional rules based on the presence of other item types in the hierarchy in which the
        target item resides. The function first finds the root of the hierarchy in which the target
        item resides, and then traverses the complete hierarchy to build a set containing the item
        types of all items present in the hierarchy, and the item types of all their compile-time
        ascendants. Finally the function calls the Schema::allowedDescendants() function to produce
        the result. */
    vector<string> allowedDescendants();
};

////////////////////////////////////////////////////////////////////

#endif

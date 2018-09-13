/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ENUMPROPERTYHANDLER_HPP
#define ENUMPROPERTYHANDLER_HPP

#include "PropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties of enumeration types. */
class EnumPropertyHandler : public PropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using PropertyHandler::PropertyHandler;

    // ================== Overriding base class functions ==================

public:
    /** Returns true if the given string can be successfully converted to a value of the property's
        type. For enumeration properties, the function returns true if the specified string matches
        one of the enumeration names for the handled property; otherwise it returns false. */
    bool isValidValue(string value) const override;

    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. For enumeration
        properties, the function inserts a name obtained by concatenating the target property's
        name and the enumeration name for the current value. In addition, the function inserts the
        names provided in the conditional expression of the "insert" attribute of the target
        property, if any. */
    virtual void insertNames() override;

    /** Accepts the specified property handler visitor. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Specific functions for this property type ==================

public:
    /** Returns a list of all enumeration names defined for the type of the handled property. */
    vector<string> values() const;

    /** Returns a list of the titles corresponding to all enumeration names defined for the type of
        the handled property, in the order corresponding to values(). */
    vector<string> titlesForValues() const;

    /** Returns the enumeration name corresponding to the default value for the handled property, or
        the empty string if unavailable. */
    string defaultValue() const;

    /** Returns the enumeration name corresponding to the value of the handled property in the
        target item. */
    string value() const;

    /** Returns the title corresponding to the value of the handled property in the
        target item. If there is no title, the empty string is returned instead. */
    string titleForValue() const;

    /** Sets the value of the handled property in the target item to the value corresponding to the
        specified enumeration name. If the specified key is invalid for this property, nothing
        happens. */
    void setValue(string value);
};

////////////////////////////////////////////////////////////////////

#endif

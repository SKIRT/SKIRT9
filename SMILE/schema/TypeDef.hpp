/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TYPEDEF_HPP
#define TYPEDEF_HPP

#include "Basics.hpp"
#include "PropertyDef.hpp"
#include <list>
class Item;

////////////////////////////////////////////////////////////////////

/** The TypeDef class represents a type definition in a SMILE schema. It is essentially a structure
    that offers setters and getters for its data members. */
class TypeDef final
{
    // ================== Constructing ==================

public:
    /** The default (and only) constructor trivially initializes the type definition to a blank
        state. All string fields are empty, there are no property definitions, and the type is
        considered to be abstract (i.e. not concrete). Use the setters to further initialize all
        relevant data members. */
    TypeDef() { }

    // ================== Copying and moving ==================

public:
    /** The copy constructor is deleted. TypeDef instances can't be copied or moved because of the
        PropertyDef instances they contain. */
    TypeDef(const TypeDef&) = delete;

    /** The assignment operator is deleted. TypeDef instances can't be copied or moved because of
        the PropertyDef instances they contain. */
    TypeDef& operator=(const TypeDef&) = delete;

    // ================== Setters ==================

public:
    /** Sets the name of the type. */
    void setName(string name) { _name = name; }

    /** Sets a description for the type that can be displayed to a user. */
    void setTitle(string title) { _title = title; }

    /** Sets the name of the immediate base type (use the empty string for the root type). */
    void setBase(string base) { _base = base; }

    /** Sets a Boolean expression using names of other types that indicates whether items of the
        type are allowed in the dataset. An empty string means always allowed. */
    void setAllowedIf(string allowedIf) { _allowedIf = allowedIf; }

    /** Sets the flag indicating that this is a concrete type rather than an abstract type (which
        is the default). */
    void setConcrete() { _concrete = true; }

    /** Sets the flag indicating that properties of sub-types of this type must be listed before
        the properties of this base type. By default (i.e. if the flag is not set), the order is
        the other way around, i.e. properties of sub-types are listed after the properties of the
        base type. */
    void setSubPropertiesFirst() { _subPropertiesFirst = true; }

    /** Type definition for a function that creates an instance of an Item subclass. */
    using Instantiator = Item* (*)();

    /** Sets a pointer to the function that creates an instance of the Item subclass described by
        this type definition. This information is set only for concrete hardcoded classes (i.e. not
        when using ghost items or for hardcoded abstract classes). */
    void setInstantiator(Instantiator instantiator) { _instantiator = instantiator; }

    /** Add a new property definition to the list for this type with the given property type, name
        and title, and returns a reference to the newly constructed PropertyDef instance so that
        its further attributes can be initialized. Properties are stored in order of addition,
        which should correspond to the order of appearance in the schema definition. Ownership of
        the PropertyDef instance is passed to the receiving type definition. */
    PropertyDef& addPropertyDef(string type, string name, string title);

    // ================== Getters ==================

public:
    /** Returns the name of the type. */
    string name() const { return _name; }

    /** Returns a description for the type that can be displayed to a user. */
    string title() const { return _title; }

    /** Returns the name of the immediate base type, or the empty string for the root type. */
    string base() const { return _base; }

    /** Returns a Boolean expression using names of other types that indicates whether items of the
        type are allowed in the dataset. An empty string means always allowed. */
    string allowedIf() const { return _allowedIf; }

    /** Returns true if this is a concrete type, and false if this is an abstract type. */
    bool concrete() const { return _concrete; }

    /** Returns true if properties of sub-types of this type must be listed before the properties
        of this base type, and false if properties of sub-types must be listed after the properties
        of the base type (the default). */
    bool subPropertiesFirst() const { return _subPropertiesFirst; }

    /** Returns a pointer to the function that creates an instance of the Item subclass described
        by this type definition. This information is available only for concrete hardcoded classes.
        The function returns a null pointer when using ghost items or for hardcoded abstract
        classes). */
    Instantiator instantiator() const { return _instantiator; }

    /** Returns (a reference to a list of) the definitions of the type's properties, in order of
        appearance in the schema. Properties of base types are NOT included in the list. */
    const std::list<PropertyDef>& propertyDefs() const { return _propertyDefs; }

    // ================== Data members ==================

private:
    string _name;           // the name of the type
    string _title;          // a description for display to a user
    string _base;           // the name of the immediate base type, or the empty string for the root type
    string _allowedIf;      // a Boolean expression using names of other types, or the empty string
    bool _concrete{false};  // true if this is a concrete type
    bool _subPropertiesFirst{false};  // true if sub-type properties must be listed first
    std::list<PropertyDef> _propertyDefs; // the definitions of the type's properties
                            // (we use std::list because it can handle items that aren't copyable)
    Instantiator _instantiator{nullptr};  // function to create Item subclass of this type
};

////////////////////////////////////////////////////////////////////

#endif

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
    TypeDef() {}

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

    /** Sets a Boolean expression that indicates whether items of the type are allowed. An empty
        string means always allowed. */
    void setAllowedIf(string allowedIf) { _allowedIf = allowedIf; }

    /** Sets a Boolean expression that indicates whether items of the type are displayed. An empty
        string means always displayed. */
    void setDisplayedIf(string displayedIf) { _displayedIf = displayedIf; }

    /** Sets a conditional value expression providing a list of extra names to be inserted in the
        global and/or local name set when an item of this type is added to the dataset, in addition
        to the names of the type and its ancestors. An empty string means that no extra names will
        be added. */
    void setInsert(string insert) { _insert = insert; }

    /** Sets the flag indicating that this is a concrete type rather than an abstract type (which
        is the default). */
    void setConcrete() { _concrete = true; }

    /** Sets the index in the property list where properties of subtypes should be listed to the
        number of properties currently held by this type definition. This has the effect of listing
        properties of subtypes just before the next property added to the type definition (or at
        the end if no further properties are added). */
    void setSubPropertyIndexHere() { _subPropertyIndex = _propertyDefs.size(); }

    /** Sets the index in the property list where properties of subtypes should be listed. If the
        index is nonnegative, properties of subtypes must be listed just before the property with
        the specified index. If the index is negative, properties of subtypes must be listed after
        the properties of the base type. This is the default behavior if this function is never
        called. */
    void setSubPropertyIndex(int index) { _subPropertyIndex = index; }

    /** Add a new property definition to the list for this type with the given property type, name
        and title, and returns a reference to the newly constructed PropertyDef instance so that
        its further attributes can be initialized. Properties are stored in order of addition,
        which should correspond to the order of appearance in the schema definition. Ownership of
        the PropertyDef instance is passed to the receiving type definition. */
    PropertyDef& addPropertyDef(string type, string name, string title);

    /** Type definition for a function that creates an instance of an Item subclass. */
    using Instantiator = Item* (*)();

    /** Sets a pointer to the function that creates an instance of the Item subclass described by
        this type definition. This information is set only for concrete hardcoded classes (i.e. not
        when using ghost items or for hardcoded abstract classes). */
    void setInstantiator(Instantiator instantiator) { _instantiator = instantiator; }

    // ================== Getters ==================

public:
    /** Returns the name of the type. */
    string name() const { return _name; }

    /** Returns a description for the type that can be displayed to a user. */
    string title() const { return _title; }

    /** Returns the name of the immediate base type, or the empty string for the root type. */
    string base() const { return _base; }

    /** Returns a Boolean expression that indicates whether items of the type are allowed. An empty
        string means always allowed. */
    string allowedIf() const { return _allowedIf; }

    /** Returns a Boolean expression that indicates whether the type is displayed. An empty string
        means always displayed. */
    string displayedIf() const { return _displayedIf; }

    /** Returns a conditional value expression providing a list of extra names to be inserted in
        the global and/or local name set when an item of this type is added to the dataset, in
        addition to the names of the type and its ancestors. An empty string means that no extra
        names will be added. */
    string insert() const { return _insert; }

    /** Returns true if this is a concrete type, and false if this is an abstract type. */
    bool concrete() const { return _concrete; }

    /** Returns the index in the property list where properties of subtypes should be listed. If
        the index is nonnegative, properties of subtypes must be listed just before the property
        with the specified index. If the index is negative, properties of subtypes must be listed
        after the properties of the base type (the default). */
    int subPropertyIndex() const { return _subPropertyIndex; }

    /** Returns the number of properties defined for the type. Properties of base types are NOT
        counted towards this number. */
    int numSubProperties() const { return _propertyDefs.size(); }

    /** Returns (a reference to a list of) the definitions of the type's properties, in order of
        appearance in the schema. Properties of base types are NOT included in the list. */
    const std::list<PropertyDef>& propertyDefs() const { return _propertyDefs; }

    /** Returns a pointer to the function that creates an instance of the Item subclass described
        by this type definition. This information is available only for concrete hardcoded classes.
        The function returns a null pointer when using ghost items or for hardcoded abstract
        classes). */
    Instantiator instantiator() const { return _instantiator; }

    // ================== Data members ==================

private:
    string _name;               // the name of the type
    string _title;              // a description for display to a user
    string _base;               // the name of the immediate base type, or the empty string for the root type
    string _allowedIf;          // a Boolean expression indicating whether items of this type are allowed
    string _displayedIf;        // a Boolean expression indicating whether this type is displayed
    string _insert;             // a conditional value expression providing a list of extra names to be inserted
    bool _concrete{false};      // true if this is a concrete type
    int _subPropertyIndex{-1};  // if nonnegative, subproperties are inserted before the property
                                // with the specified index rather than added at the end of the list
    std::list<PropertyDef> _propertyDefs;  // the definitions of the type's properties
                                           // (we use std::list because it can handle items that aren't copyable)
    Instantiator _instantiator{nullptr};   // function to create Item subclass of this type
};

////////////////////////////////////////////////////////////////////

#endif

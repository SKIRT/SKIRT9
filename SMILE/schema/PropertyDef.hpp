/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROPERTYDEF_HPP
#define PROPERTYDEF_HPP

#include "Basics.hpp"
class PropertyAccessor;

////////////////////////////////////////////////////////////////////

/** The PropertyDef class represents SMILE schema property definitions. It is essentially a
    structure that offers setters and getters for its data members. A PropertyDef instance can hold
    the definition for any property type; attributes that are not relevant for a particular
    property type are simply ignored by the corresponding property handler. */
class PropertyDef final
{
    // ================== Constructing ==================

public:
    /** The constructor sets the property definition's type (one of the supported property type
        names, e.g. "IntProperty", "ItemListProperty"), its name (for use by a program) and its
        description (for display to a user). All other property attributes are initialized to the
        value that corresponds to an empty string or container (depending on its data type). Use
        the setters to further initialize any relevant data members as desired. */
    PropertyDef(string type, string name, string title);

    /** The destructor releases any resources being held, such as the accessor info block. */
    ~PropertyDef();

    // ================== Copying and moving ==================

public:
    /** The copy constructor is deleted. PropertyDef instances can't be copied or moved because
        they contain a unique pointer to a PropertyAccessor instance. */
    PropertyDef(const PropertyDef&) = delete;

    /** The assignment operator is deleted. PropertyDef instances can't be copied or moved because
        they contain a unique pointer to a PropertyAccessor instance. */
    PropertyDef& operator=(const PropertyDef&) = delete;

    // ================== Setters ==================

public:
    /** Sets the Boolean expression indicating whether this property is relevant. If the Boolean
        expression is empty, the property is always relevant. */
    void setRelevantIf(string value) { _relevantIf = value; }

    /** Sets the Boolean expression indicating whether this property is displayed. If the Boolean
        expression is empty, the property is always displayed. */
    void setDisplayedIf(string value) { _displayedIf = value; }

    /** Sets the Boolean expression indicating whether the property is required. If the Boolean
        expression is empty, the property is always required. */
    void setRequiredIf(string value) { _requiredIf = value; }

    /** Sets a conditional value expression providing a list of extra names to be inserted in the
        global and/or local name set when a value of this item is entered into the dataset, in
        addition to the name automatically associated with the property. An empty string means that
        no extra names will be added. */
    void setInsert(string insert) { _insert = insert; }

    /** Sets the conditional value expression defining the default value for the property in case
        the property is missing. If the string is empty, there is no default value. */
    void setDefaultValue(string value) { _default = value; }

    /** Sets the minimum value for the property, in a format compatible with the particular
        property type. If the string is empty, the minimum value is determined by the
        implementation of the property type. */
    void setMinValue(string value) { _min = value; }

    /** Sets the maximum value for the property, in a format compatible with the particular
        property type. If the string is empty, the maximum value is determined by the
        implementation of the property type. */
    void setMaxValue(string value) { _max = value; }

    /** Sets the name of the physical quantity represented by this property. This must match one of
        the quantity names defined in this schema and determines the allowed and default units for
        the property values. If the string is empty, the property values are dimensionless. */
    void setQuantity(string value) { _quantity = value; }

    /** Sets the name of the base type for the instances of this property. This must match one of
        the type names defined in this schema. This information is relevant only for compound
        properties (i.e. properties containing other items of some type), and for those properties
        it should not be left empty. */
    void setBase(string value) { _base = value; }

    /** Adds a name/title pair for one of the enumeration values defined for this property.
        Enumeration values must be added in order of occurrence in the schema definition. This
        information is relevant only for enumeration properties, and for those properties the list
        of enumeration values should not be left empty. */
    void addEnumeration(string name, string title);

    /** Hands a property accessor block to the property definition; ownership is passed to the
        property definition. */
    void setAccessor(const PropertyAccessor* accessor);

    // ================== Getters ==================

public:
    /** Returns the type of the property as one of the supported property type names, e.g.
        "IntProperty", "ItemListProperty". */
    string type() const { return _type; }

    /** Returns the name of the property for use by a program. */
    string name() const { return _name; }

    /** Returns the description for the property for display to a user. */
    string title() const { return _title; }

    /** Returns the Boolean expression indicating whether this property is relevant. If the Boolean
        expression is empty, the property is always relevant. */
    string relevantIf() const { return _relevantIf; }

    /** Returns the Boolean expression indicating whether this property is displayed. If the
        Boolean expression is empty, the property is always displayed. */
    string displayedIf() const { return _displayedIf; }

    /** Returns the Boolean expression indicating whether this property is required. If the
        Boolean expression is empty, the property is always required. */
    string requiredIf() const { return _requiredIf; }

    /** Returns a conditional value expression providing a list of extra names to be inserted in the
        global and/or local name set when a value of this item is entered into the dataset, in
        addition to the name automatically associated with the property. An empty string means that
        no extra names will be added. */
    string insert() const { return _insert; }

    /** Returns the default value for the property in case the property is missing, in a format
        compatible with the particular property type. If the string is empty, there is no default
        value. */
    string defaultValue() const { return _default; }

    /** Returns the minimum value for the property, in a format compatible with the particular
        property type. If the string is empty, the minimum value is determined by the
        implementation of the property type. */
    string minValue() const { return _min; }

    /** Returns the maximum value for the property, in a format compatible with the particular
        property type. If the string is empty, the maximum value is determined by the
        implementation of the property type. */
    string maxValue() const { return _max; }

    /** Returns the name of the physical quantity represented by this property. This must match one
        of the quantity names defined in this schema and determines the allowed and default units
        for the property values. If the string is empty, the property values are dimensionless. */
    string quantity() const { return _quantity; }

    /** Returns the name of the base type for the instances of this property. This must match one
        of the type names defined in this schema. This information is relevant only for compound
        properties (i.e. properties containing other items of some type). */
    string base() const { return _base; }

    /** Returns a list of all enumeration names defined for this property, in order of occurrence
        in the schema. For properties other than enumeration properties, this function returns an
        empty list. */
    const vector<string>& enumNames() const { return _enumNames; }

    /** Returns a list of the titles corresponding to all enumeration names defined for this
        property, in the order corresponding to the list returned by enumNames(). For properties
        other than enumeration properties, this function returns an empty list. */
    const vector<string>& enumTitles() const { return _enumTitles; }

    /** Returns the title corresponding to the specified enumeration name, as defined for this
        property. If the specified enumeration name is not defined for the property, the function
        returns the empty string. */
    string enumTitle(string enumName) const;

    /** Returns a pointer to the property accessor block for the property definition, or nullptr if
        no accessor block was set; ownership stays with the propery definition. */
    const PropertyAccessor* accessor() const { return _accessor.get(); }

    // ================== Data members ==================

private:
    string _type;                // property type
    string _name;                // name of the property
    string _title;               // description of the property
    string _relevantIf;          // Boolean expression indicating whether this property is relevant
    string _displayedIf;         // Boolean expression indicating whether the property should be displayed
    string _requiredIf;          // Boolean expression indicating whether the property is required
    string _insert;              // conditional value expression providing a list of extra names to be inserted
    string _default;             // default value of the property (as a conditional value expression)
    string _min;                 // minimum value of the property
    string _max;                 // maximum value of the property
    string _quantity;            // name of the physical quantity represented by this property
    string _base;                // name of the base type for the instances of this property
    vector<string> _enumNames;   // list of names for the enumeration values of this property
    vector<string> _enumTitles;  // list of titles for the enumeration values of this property, in the same order
    std::unique_ptr<const PropertyAccessor> _accessor;  // accessor block for the property's getters and setters
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROPERTYHANDLER_HPP
#define PROPERTYHANDLER_HPP

#include "Basics.hpp"
class Item;
class NameManager;
class PropertyDef;
class PropertyHandlerVisitor;
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** PropertyHandler is the abstract base class for handling a property of a SMILE data item. There
    is a specific subclass for each supported property type. A PropertyHandler instance combines
    knowledge about the schema definition describing the SMILE dataset being handled, the target
    data item within that dataset, the target property within that item, and the type of that
    property. The latter aspect is baked into the code of the PropertyHandler subclass in use,
    while the first three aspects are provided during construction. */
class PropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. The last argument, a pointer to the name manager used for the current
        dataset, allows the property handler to access the global and local name sets for
        evaluating Boolean expressions.

        The caller must guarantee that the specified instances are not deleted during the lifetime
        of the property handler (the property handler does not assume ownership, nor is there a
        reference counting scheme). */
    PropertyHandler(Item* target, const PropertyDef* property, const SchemaDef* schema, NameManager* nameMgr);

    /** Destructs the property handler. This virtual destructor is declared here because
        PropertyHandler is the top-level class in the hierarchy of property handlers. */
    virtual ~PropertyHandler() = default;

    /** The copy constructor is deleted because instances of derived classes should never be copied
        or moved. */
    PropertyHandler(const PropertyHandler&) = delete;

    /** The assignment operator is deleted because instances of derived classes should never be
        copied or moved. */
    PropertyHandler& operator=(const PropertyHandler&) = delete;

    // ================== Getters ==================

protected:
    /** Returns the SMILE data item for which a property is being handled. */
    Item* target() const { return _target; }

    /** Returns the property definition for the SMILE property being handled. */
    const PropertyDef* property() const { return _property; }

public:
    /** Returns the schema definition describing the SMILE dataset in which the target item
        resides. */
    const SchemaDef* schema() const { return _schema; }

    /** Returns the name manager for the SMILE dataset in which the target item resides. */
    NameManager* nameManager() const { return _nameMgr; }

    /** Returns the type of the target item. */
    string type() const;

    /** Returns the name of the handled property. */
    string name() const;

    /** Returns the title (used for display to a user) for the handled property. */
    string title() const;

    /** Returns true if the handled property is relevant in the current dataset configuration. A
        property that does not have the "relevantIf" attribute is always relevant. A property that
        has the "relevantIf" attribute is relevant if the provided Boolean expression evaluates to
        true, after replacing names currently present in the global or local name sets by true, and
        other names by false. */
    bool isRelevant() const;

    /** Returns true if the handled property should be displayed in the current dataset
        configuration. A property that does not have the "displayedIf" attribute is always
        displayed. A property that has the "displayedIf" attribute is displayed if the provided
        Boolean expression evaluates to true, after replacing names currently present in the global
        or local name sets by true, and other names by false. */
    bool isDisplayed() const;

    /** Returns true if the handled property is required to have a nonempty value in the current
        dataset configuration. A property that does not have the "requiredIf" attribute is always
        required. A property that has the "requiredIf" attribute is required if the provided
        Boolean expression evaluates to true, after replacing names currently present in the global
        or local name sets by true, and other names by false. The "requiredIf" attribute (and thus
        the return value of this function) is intended for complex property types including
        strings, compound types, and lists. Scalar numeric and enumeration property types ignore
        its value. */
    bool isRequired() const;

    /** Returns true if the handled property has a valid default value in the current dataset
        configuration, or false if not. A property that does not have the "default" attribute does
        not have a default value. A property that has the "default" attribute has a default value
        if the provided conditional value expression evaluates to a valid value for the property
        type, after replacing names currently present in the global or local name sets by true, and
        other names by false. */
    bool hasDefaultValue() const;

    /** Returns true if the given string can be successfully converted to a value of the property's
        type. For an empty string, the function always returns false. This function must be
        overridden by property handler subclasses. */
    virtual bool isValidValue(string value) const = 0;

    /** Returns true if the handled property type is compound in the sense that the property may
        hold (in other words, aggregate) other items that are part of the item hierarchy. The
        default implementation in this abstract class returns false, which is appropriate for
        plain, non-compound property types. This function must be overridden by property handler
        subclasses that implement compound property types. */
    virtual bool isCompound() const;

    /** Returns true if the value of the target property has been modified by this handler;
        otherwise returns false. */
    bool hasChanged() const;

    // ================== Setters & Modifiers ==================

protected:
    /** Sets the flag indicating that the value of the target property has been modified by this
        handler, and invokes the insertNames() function. Subclasses must call this function when
        the target property is updated. */
    void setChanged();

public:
    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. Because the list
        of inserted names depends on the property type, this function must be implemented in every
        subclass. */
    virtual void insertNames() = 0;

    /** Accepts the specified property handler visitor. This function is part of the "visitor"
        design pattern implementation used to handle properties of various types. It must be
        implemented in every subclass. */
    virtual void acceptVisitor(PropertyHandlerVisitor* visitor) = 0;

    // ================== Utilities for editing ==================

public:
    /** This function stores a flag in the target dataset item to indicate whether the user has
        configured the property being handled during this session, depending on the value of the
        specified argument. If the argument is omitted, the default value of true is used. */
    void setConfigured(bool configured=true);

    /** This function returns true if the setPropertyConfigured() function was last called with a
        true argument during this session for this combination of target dataset item and property,
        and false otherwise. In other words, it returns true if the property being handled has been
        configured by the user for the target dataset item. */
    bool isConfigured();

    // ================== Data members ==================

private:
    Item* _target{nullptr};                 // the item being handled
    const PropertyDef* _property{nullptr};  // the property definition for the property being handled
    const SchemaDef* _schema{nullptr};      // the schema definition for the dataset in which the item resides
    NameManager* _nameMgr{nullptr};         // the name manager for the dataset in which the item resides
    bool _changed{false};           // becomes true if the target property has been modified by this handler
};

////////////////////////////////////////////////////////////////////

#endif

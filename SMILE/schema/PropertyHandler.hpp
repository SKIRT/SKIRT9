/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROPERTYHANDLER_HPP
#define PROPERTYHANDLER_HPP

#include "Basics.hpp"
class Item;
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
        schema definition. The caller must guarantee that the specified instances are not deleted
        during the lifetime of the property handler (the property handler does not assume
        ownership, nor is there a reference counting scheme). */
    PropertyHandler(Item* target, const PropertyDef* property, const SchemaDef* schema);

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

    /** Returns the property definition for the property being handled. */
    const PropertyDef* property() const { return _property; }

public:
    /** Returns the schema definition describing the SMILE dataset in which the target item
        resides. */
    const SchemaDef* schema() const { return _schema; }

    /** Returns the type of the target item. */
    string type() const;

    /** Returns the name of the handled property. */
    string name() const;

    /** Returns the title (used for display to a user) for the handled property. */
    string title() const;

    /** Returns true if the handled property is silent, i.e. an interactive Q&A session should
        not ask for a value corresponding to this property. A silent property must have an
        appropriate default value. */
    bool isSilent() const;

    /** Returns true if the handled property has a non-empty "RelevantIf" attribute so that it may
        become irrelevant depending on the values of other properties. */
    bool hasRelevantIf() const;

    /** Returns true if the handled property is relevant given the current values of the other
        properties in the target item, and false otherwise. A property that does not have the
        "RelevantIf" attribute is always relevant. A property that has the "RelevantIf" attribute
        is relevant if the provided Boolean expression evaluates to true. For this evaluation, a
        property name in the expression are replaced by true if that property is relevant in its
        own right (as determined by a recursive call to isRelevant()), and has a value that would
        test as true in a condition (as determined by isTrueInCondition()). Otherwise the property
        name is replaced by false in the expression. */
    bool isRelevant() const;

    /** Returns true if the value of the handled property would test as true in a condition, and
        false if it would test as false or if the property type can't be used in a condition. For
        example a boolean property would simply return its value; an integer property would return
        true if its value is nonzero. This function is used by isRelevant() to test a simple
        condition on the value of a property in the current dataset in a way that doesn't depend
        on the property type. The default implementation provided in this abstract class always
        returns false; this function must be overridden by property handler subclasses that support
        condition testing. */
    virtual bool isTrueInCondition() const;

    /** Returns true if the handled property is optional (i.e. its value may remain unset), or
        false if not. The default implementation in this abstract class returns false; this
        function must be overridden by property handler subclasses that support optional values. */
    virtual bool isOptional() const;

    /** Returns true if the handled property has a valid default value, or false if not. The
        default implementation in this abstract class returns false; this function must be
        overridden by property handler subclasses that may offer a default value. */
    virtual bool hasDefaultValue() const;

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
        handler. Subclasses must call this function when the target property is updated. */
    void setChanged();

public:
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
    bool _changed{false};           // becomes true if the target property has been modified by this handler
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITEMPROPERTYHANDLER_HPP
#define ITEMPROPERTYHANDLER_HPP

#include "AbstractItemPropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties of type "pointer to item". */
class ItemPropertyHandler : public AbstractItemPropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using AbstractItemPropertyHandler::AbstractItemPropertyHandler;

    // ================== Overriding base class functions ==================

public:
    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. For item
        properties, the function inserts the target property's name if the current property value
        is not null (i.e. an item is present), and does not insert any names if the value is null
        (i.e no item is present). */
    void insertNames() override;

    /** Accepts the specified visitor. This function is part of the "visitor" design pattern
        implementation used to handle properties of various types. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Specific functions for this property type ==================

public:
    /** Returns the value of the handled property in the target item. There is no transfer of
        ownership. */
    Item* value() const;

    /** Sets the value of the handled property in the target item so that it points to the
        specified item. The target item assumes ownership of the specified item instance. The
        function returns false if the property couldn't be set (e.g. because the specified item has
        an inappropriate type). */
    bool setValue(Item* value);

    /** Constructs a new item instance of the specified type and sets the value of the handled
        property in the target item so that it points to this new instance. The target item assumes
        ownership of the new instance. The function returns false if the property couldn't be set
        (e.g. because the specified item type is inappropriate). */
    bool setToNewItemOfType(string type);

    /** Sets the value of the handled property in the target item to a null pointer, removing any
        previously owned item instance. */
    void setToNull();
};

////////////////////////////////////////////////////////////////////

#endif

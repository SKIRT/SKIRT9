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
    /** Returns a list of SMILE data items including the immediate children of the target item in
        the dataset in which the target item resides. The implementation in this class returns the
        result of the value() function encapsulated in a list, or the empty list if the value is
        null. */
    vector<Item*> children() const override;

    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. For item
        properties, the function inserts the target property's name if the current property value
        is not null (i.e. an item is present), and does not insert any names if the value is null
        (i.e no item is present). In addition, the function inserts the names provided in the
        conditional expression of the "insert" attribute of the target property, if any.

        Furthermore, if the current property value is not null (i.e. an item is present), the
        function inserts the name of the item's type and of all its ancestor types, recursively. In
        addition, it inserts the names provided in the conditional expression of the "insert"
        attribute of these types, if any. */
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

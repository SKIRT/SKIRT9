/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITEMLISTPROPERTYHANDLER_HPP
#define ITEMLISTPROPERTYHANDLER_HPP

#include "AbstractItemPropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** This class handles SMILE data item properties of type "list of pointer to item". */
class ItemListPropertyHandler : public AbstractItemPropertyHandler
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
        result of the value() function. */
    vector<Item*> children() const override;

    /** Causes the name manager associated with this handler to insert names into the global and/or
        local name sets corresponding to the current value of the target property. For item list
        properties, the function inserts the target property's name if the current property value
        is a nonempty list (i.e. including at least one item), and does not insert any names if the
        value is the empty list. In addition, the function inserts the names provided in the
        conditional expression of the "insert" attribute of the target property, if any.

        Furthermore, for each of the items in the current property value, the function inserts the
        name of the item's type and of all its ancestor types, recursively. In addition, it inserts
        the names provided in the conditional expression of the "insert" attribute of these types,
        if any. */
    void insertNames() override;

    /** Accepts the specified visitor. This function is part of the "visitor" design pattern
        implementation used to handle properties of various types. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Specific functions for this property type ==================

public:
    /** Returns the value of the handled property in the target item. */
    vector<Item*> value() const;

    /** Empties the list held by the handled property in the target item, removing any previously
        owned item instances. */
    void setToEmpty();

    /** Adds the specified item to the list held by the handled property in the target
        item. The target item assumes ownership of the specified instance. The function returns
        false if the item couldn't be added (e.g. because it has an inappropriate type). */
    bool addValue(Item* value);

    /** Constructs a new item instance of the specified type and adds it to the list held by the
        handled property in the target item. The target item assumes ownership of the new instance.
        The function returns false if a new item couldn't be added (e.g. because the specified type
        is inappropriate). */
    bool addNewItemOfType(string type);

    /** Inserts the specified item at the specified index into the list held by the handled
        property in the target item. The target item assumes ownership of the specified instance.
        The function returns false if the item couldn't be added (e.g. because it has an
        inappropriate type). */
    bool insertValue(int index, Item* value);

    /** Constructs a new instance of the specified type and inserts it at the specified index into
        the list held by the handled property in the target item. The target item assumes ownership
        of the new instance. The function returns false if a new item couldn't be added (e.g.
        because the specified item type is inappropriate). */
    bool insertNewItemOfType(int index, string type);

    /** Removes the item with the specified zero-based index from the list held by the handled
        property in the target item. The removed item is deleted. The function returns false if the
        item couldn't be removed (e.g. because the index is out of range). */
    bool removeValueAt(int index);

    // ================== Utilities for editing ==================

public:
    /** This function stores the selected row index for the handled property in the target item to
        the specified integer value. */
    void storeSelectedRow(int row);

    /** This function returns the stored selected row index for the handled property in the target
        item. If the selected row index has never been stored for this property and item, the
        function returns zero. */
    int retrieveSelectedRow();
};

////////////////////////////////////////////////////////////////////

#endif

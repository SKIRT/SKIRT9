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
    /** Returns true if the list held by the handled property is not empty, and false if it is empty. */
    bool isTrueInCondition() const override;

    /** Accepts the specified visitor. This function is part of the "visitor" design pattern
        implementation used to handle properties of various types. */
    void acceptVisitor(PropertyHandlerVisitor* visitor) override;

    // ================== Functionality for this property type ==================

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

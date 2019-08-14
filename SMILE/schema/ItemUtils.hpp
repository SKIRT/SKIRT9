/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITEMUTILS_HPP
#define ITEMUTILS_HPP

#include "Basics.hpp"
class Item;
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** This class offers a set of static utility functions for working with SMILE dataset items in the
    context of a wizard that needs to move back and forth in the hierarchy. */
class ItemUtils final
{
public:
    /** This function sets the \em configured state for the specified property in the specified
        item to the specified integer value. A newly created item has a \em configured state of
        zero. */
    static void setPropertyConfiguredState(Item* item, string property, int configured);

    /** This function sets the \em configured state of all properties in the specified dataset
        hierarchy to one. The function calls itself recursively to process the children of the
        specified root item. */
    static void setHierarchyConfigured(const SchemaDef* schema, Item* root);

    /** This function returns the \em configured state for the specified property in the specified
        item. If the \em configured state has never been set for this property and item, the
        function returns zero. */
    static int propertyConfiguredState(Item* item, string property);

    /** This function sets the \em complete state for the specified item to true. A newly created
        item has a \em complete state of false. */
    static void setItemComplete(Item* item);

    /** This function sets the \em complete state for all items in the specified dataset hierarchy
        to true. The function calls itself recursively to process the children of the specified
        root item. */
    static void setHierarchyComplete(Item* root);

    /** This function clears the \em complete state for the specified item and for all its
        ascendants in the run-time hierarchy. */
    static void setItemIncomplete(Item* item);

    /** This function returns the \em complete state for the specified item. If the \em complete
        state has never been set for this item, the function returns false. */
    static bool isItemComplete(Item* item);

    /** This function stores the selected row index for the specified property in the specified
        item to the specified integer value. The function should be called only for item list
        properties, but the current implementation does not enforce this. */
    static void storeSelectedRow(Item* item, string property, int row);

    /** This function returns the stored selected row index for the specified property in the
        specified item. If the selected row index has never been stored for this property and item,
        the function returns zero. */
    static int retrieveSelectedRow(Item* item, string property);
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GHOSTITEM_HPP
#define GHOSTITEM_HPP

#include "Item.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

/** An instance of the GhostItem class represents an arbitrary item in any dataset described by a
    SMILE schema. The GhostItem class does not provide any actual item-specific functionality: the
    intention is to represent and manipulate a SMILE dataset in an abstract way without expecting
    it to actually do something. To this end, the GhostItem class implements the setters and
    getters for property values of the various supported property types, just like the Item base
    class. However, the property values are stored in generic dictionaries rather than in custom
    data members tailored to the type being represented. */
class GhostItem final : public Item
{
    // ================== Constructing ==================

public:
    /** The constructor sets the item's type, which should match one of the concrete types defined
        in the SMILE schema describing the dataset to which the item belongs. The item's parent is
        initialized to the null pointer. For all items other than the root node, the parent is set
        later on by handing the item to an item or item list property of another item (which calls
        the private setParent() function). The constructor further initializes the item's state so
        that no properties have been setup, which also implies that the item has no children. */
    GhostItem(string type) : _type(type) {}

    // ================== Property Setters ==================

public:
    /** Sets the specified string property to the specified value. */
    void setStringProperty(const PropertyDef* property, string value) override;

    /** Sets the specified Boolean property to the specified value. */
    void setBoolProperty(const PropertyDef* property, bool value) override;

    /** Sets the specified integer property to the specified value. */
    void setIntProperty(const PropertyDef* property, int value) override;

    /** Sets the specified enumeration property to the specified value. */
    void setEnumProperty(const PropertyDef* property, string value) override;

    /** Sets the specified floating point property to the specified value. */
    void setDoubleProperty(const PropertyDef* property, double value) override;

    /** Sets the specified floating point list property to the specified value. */
    void setDoubleListProperty(const PropertyDef* property, vector<double> value) override;

    /** Sets the specified item property to the specified item, and sets the parent of the
        specified item to the receiving item so that ownership is handed over to the receiving
        item. If the specified property already holds another item, that item is destructed. */
    void setItemProperty(const PropertyDef* property, Item* item) override;

    /** Sets the specified item list property to the empty list, removing any items that were present. */
    void clearItemListProperty(const PropertyDef* property) override;

    /** Adds the specified item to the specified item list property, and sets the parent of the
        specified item to the receiving item so that ownership is handed over to the receiving
        item. Assume that the item list property has N items when this function is called (i.e.
        before anything is changed). If the specified zero-based index is in the range [0,N-1], the
        new item is inserted at the specified index (i.e. just before the item that was there
        previously). If the index is out of this range (i.e. it is negative or at least N), the new
        item is added at the end of the list. */
    void insertIntoItemListProperty(const PropertyDef* property, int index, Item* item) override;

    /** Removes the item at the specified index from the specified item list property and deletes
        it. If the index is out of range, nothing happens. */
    void removeFromItemListProperty(const PropertyDef* property, int index) override;

    // ================== Property Getters ==================

public:
    /** Returns the name of the item's type. */
    string type() const override;

    /** Returns the value of the specified string property. If the property does not exist, the
        function throws an error. */
    string getStringProperty(const PropertyDef* property) const override;

    /** Returns the value of the specified Boolean property. If the property does not exist, the
        function throws an error. */
    bool getBoolProperty(const PropertyDef* property) const override;

    /** Returns the value of the specified integer property. If the property does not exist, the
        function throws an error. */
    int getIntProperty(const PropertyDef* property) const override;

    /** Returns the value of the specified enumeration property. If the property does not exist,
        the function throws an error. */
    string getEnumProperty(const PropertyDef* property) const override;

    /** Returns the value of the specified floating point property. If the property does not exist,
        the function throws an error. */
    double getDoubleProperty(const PropertyDef* property) const override;

    /** Returns the value of the specified floating point list property. If the property does not
        exist, the function throws an error. */
    vector<double> getDoubleListProperty(const PropertyDef* property) const override;

    /** Returns the value of the specified item property. Ownership for the returned item is \em
        not handed over. If the property does not exist, the function throws an error. */
    Item* getItemProperty(const PropertyDef* property) const override;

    /** Returns the value of the specified item list property. Ownership for the returned items is
        \em not handed over. If the property does not exist, the function throws an error. */
    vector<Item*> getItemListProperty(const PropertyDef* property) const override;

    // ================== Data members ==================

private:
    string _type;  // item type

    // containers for properties of each supported type:  <property-name, property-value>
    std::unordered_map<string, string> _stringProperties;
    std::unordered_map<string, bool> _boolProperties;
    std::unordered_map<string, int> _intProperties;
    std::unordered_map<string, string> _enumProperties;
    std::unordered_map<string, double> _doubleProperties;
    std::unordered_map<string, vector<double>> _doubleListProperties;
    std::unordered_map<string, Item*> _itemProperties;
    std::unordered_map<string, vector<Item*>> _itemListProperties;
};

////////////////////////////////////////////////////////////////////

#endif

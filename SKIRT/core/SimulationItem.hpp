/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATIONITEM_HPP
#define SIMULATIONITEM_HPP

#include "ItemInfo.hpp"
#include <typeinfo>

////////////////////////////////////////////////////////////////////

/** This is the base class for all classes in the simulation hierarchy. It inherits from Item so
    that we can use its facilities for object hierarchies and for class introspection. Simulation
    items form a compile-time class hierarchy through inheritance (with the SimulationItem class at
    the top), and a run-time object hierarchy using the Item parent/children mechanism (usually
    with an instance of a MonteCarloSimulation subclass at the top). */
class SimulationItem : public Item
{
    ITEM_ABSTRACT(SimulationItem, Item, "a simulation item")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function performs setup for this simulation item and for all its descendents.
        Specifically it invokes setupSelfBefore() on itself, then setup() on all its children, and
        finally setupSelfAfter() on itself. As a rule, the constructor of a SimulationItem subclass
        sets any attributes to null values rather than to meaningful defaults (to avoid the extra
        construction time). Thus all attributes in the simulation hierarchy must have been
        explicitly set by the caller before invoking setup(). If setup() has already been invoked
        for the same item, this function does nothing. Do not override this function in subclasses;
        implement setupSelfBefore() and/or setupSelfAfter() instead. */
    void setup();

protected:
    /** This function should be overridden in a subclass to verify validity and completeness of the
        attribute values set for the class instance, and to perform any setup that must happen \em
        before the children of the class instance have been setup. This implementation here in the
        SimulationItem class does nothing. Any implementation overriding this function should start
        by calling the same function in its immediate base class. */
    virtual void setupSelfBefore();

    /** This function should be overridden in a subclass to perform any setup that must happen \em
        after the children of the class instance have been setup. This implementation here in the
        SimulationItem class does nothing. Any implementation overriding this function should start
        by calling the same function in its immediate base class. */
    virtual void setupSelfAfter();

    //======================== Other Functions =======================

public:
    /** This function returns the concatenation the type of the simulation item and its
        human-readable name if it has one. The resulting string can be used, for example, in log
        messages to identify the item and differentiate it from other items of a similar type. */
    string typeAndName() const;

    /** This function returns a human-readable name for the simulation item, or the empty string if
        it does not have one. The item name can be used to differentiate the item from other items
        of a similar type. The implementation in this base class returns the empty string.
        Subclasses can override the implementation to return something meaningful. */
    virtual string itemName() const;

    //======================== Finding items =======================

public:
    /** This template function looks for a simulation item of the type specified as template
        argument in the hierarchy containing the receiving object and, if found, returns a pointer
        to that item. The class specified as template parameter must inherit from SimulationItem.
        The function searches for an appropriate item among all ancestors of the receiving item
        (including the receiving item itself) and their immediate children. In other words, it
        recursively runs upwards along the ancestors and goes just one level down for each
        ancestor. The function returns the first appropriate item found; if multiple items of the
        same type exist in the hierarchy, it is undefined which one of these will be returned.

        If the \em setup flag is true (the default value), the function invokes the setup()
        function on the item before it is returned; if no appropriate item is found, a FatalError
        is thrown. If the \em setup flag is false, the function does not perform setup on
        the item, and if no appropriate item is found, the function returns a null pointer. */
    template<class T> T* find(bool setup = true) const
    {
        static_assert(std::is_base_of<SimulationItem, T>::value,
                      "Requested type in find() does not inherit SimulationItem");
        return dynamic_cast<T*>(find(setup, [](Item* item) -> SimulationItem* { return dynamic_cast<T*>(item); }));
    }

    /** This template function looks for a simulation item that offers the interface specified as
        template argument in the hierarchy containing the receiving object and, if found, returns a
        pointer to that item (or rather, a pointer to the item after it has been dynamically cast
        to the requested interface type). The function returns the first appropriate item found.
        The portion of the hierarchy searched by the function for an appropriate item is defined by
        the value of the \em levels argument as follows.

        For levels=0, only the receiving object is considered. For levels=-N (with N a positive
        integer), the receiving object plus N of its ancestors are considered in the order from
        nearest to more distant ancestor. For levels=N (with N a positive integer), the receiving
        object plus N levels of children are considered. Children at the same level are searched in
        an undefined order, except that children part of the same item list property are guaranteed
        to be searched in configuration order. Multiple levels are searched depth-first, i.e. the
        children of an item are searched before the next item on the same level is considered. The
        default value of levels is -999999, which has the effect of searching the receiving object
        plus all of its ancestors.

        For a simulation item to be considered as offering the requested interface, two conditions
        must be fullfilled. First, obviously, the item's class must inherit from and implement the
        interface. Second, the offersInterface() function, when invoked on the item under
        consideration for the requested interface type, must return \em true. The default
        implementation of the offersInterface() function provided in this base class always returns
        true. Overriding it in a subclass allows the subclass to decide at run time whether the
        conditions for offering a certain interface are fullfilled.

        If the \em setup flag is true (the default value), the function invokes the setup()
        function on the item before it is returned; if no appropriate item is found, a FatalError
        is thrown. If the \em setup flag is false, the function does not perform setup on the item,
        and if no appropriate item is found, the function returns a null pointer. */
    template<class T> T* interface(int levels = -999999, bool setup = true) const
    {
        return dynamic_cast<T*>(interface(levels, setup, [](SimulationItem* item) {
            return dynamic_cast<T*>(item) != nullptr && item->offersInterface(typeid(T));
        }));
    }

private:
    /** This is the private implementation used by the find() template function. The first argument
        has the same semantics as the \em setup argument of the template function. The second
        argument accepts a function that dynamically casts a given Item instance to the requested
        simulation item type (or returns a null pointer if the Item instance is not of that type).
        */
    Item* find(bool setup, SimulationItem* castToRequestedType(Item*)) const;

    /** This is the private implementation used by the interface() template function. The first two
        arguments have the same semantics as the corresponding arguments of the template function.
        The last argument accepts a function that returns true if the given simulation item
        implements the requested interface, and false otherwise. */
    SimulationItem* interface(int levels, bool setup, bool offersRequestedInterface(SimulationItem*)) const;

protected:
    /** This function is for use only by the interface() function. After detecting that the
        receiving item implements the specified interface (i.e. its type inherits the interface),
        the interface() function invokes the offersInterface() function to ensure that the object
        can actually offer the interface.

        Thus, this function must return true if the receiving item actually offers (not just
        implements) the specified interface in the current run-time environment, and false if not.
        The default implementation of the function provided in this base class always returns true.
        Overriding it in a subclass allows the subclass to decide at run time whether the
        conditions for offering a certain interface are fullfilled. */
    virtual bool offersInterface(const std::type_info& interfaceTypeInfo) const;

    //======================== Data Members ========================

private:
    bool _setupStarted{false};  // becomes true when setup for the item begins
};

////////////////////////////////////////////////////////////////////

#endif

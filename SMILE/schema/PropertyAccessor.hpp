/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROPERTYACCESSOR_HPP
#define PROPERTYACCESSOR_HPP

#include "Basics.hpp"
class Item;

////////////////////////////////////////////////////////////////////

/** The PropertyAccessor template class hierarchy implements the heart of the technique that allows
    type-safe access to discoverable properties of Item subclasses based on property name. There
    are three inheritance levels in the property accessor class hierarchy. At level zero, the
    PropertyAccessor class serves as a generic base class for all classes in the hierarchy. It has
    a virtual destructor to enable dynamic casting for PropertyAccessor subclass instances.
    Otherwise, it does not offer any functionality at all. PropertyAccessor is not a template
    class, but all of its subclasses are.

    At inheritance level one, the TypedPropertyAccessor template class defines the pure virtual
    interface for getting and setting a property of value type V, given as template parameter. The
    template is specialized for a value type representing a list of items, because this type
    requires a different set of setters to manage individual items in the list.

    At inheritance level two, the TargetTypedPropertyAccessor template class actually implements
    the property access interface for a value type V and for a target type T (which must inherit
    Item), both given as template parameters. TargetTypedPropertyAccessor inherits
    TypedPropertyAccessor with an appropriate value type. The template is specialized for several
    value types to allow special processing, such as casting pointers between Item and one of its
    subclasses, or translating an enumeration value to an integer and vice versa.

    This setup enables a usage pattern as follows. The PROPERTY macro's in the Item.hpp header
    expand to functions that load appropriate metadata for the discoverable properties defined in
    Item subclasses. These functions are invoked at program startup by the item registry. Each of
    these functions creates a TargetTypedPropertyAccessor instance with the value type of the
    property and the target type of the Item subclass in which the property is defined. The item
    registry subsequently stores this property accessor instance in a PropertyDef instance using
    its generic PropertyAccessor type. The generic property access functions defined in the Item
    class, e.g. getIntProperty(), retrieve the PropertyAccessor instance from the property
    definition, cast it to a TypedPropertyAccessor instance of the appropriate value type, and
    invoke the latter's access functions to actually access the property in the target item. */
class PropertyAccessor
{
protected:
    PropertyAccessor() {}

public:
    PropertyAccessor(const PropertyAccessor&) = delete;
    PropertyAccessor& operator=(const PropertyAccessor&) = delete;
    virtual ~PropertyAccessor() {}
};

////////////////////////////////////////////////////////////////////

/** The TypedPropertyAccessor template class sits at level one of the PropertyAccessor class
    hierarchy. It defines the pure virtual interface for getting and setting a property of value
    type V, given as template parameter. The template is specialized for a value type representing
    a list of items, because this type requires a different set of setters to manage individual
    items in the list. Refer to the PropertyAccessor class for more information. */
template<typename V> class TypedPropertyAccessor : public PropertyAccessor
{
protected:
    TypedPropertyAccessor() {}

public:
    virtual V targetValue(const Item* target) const = 0;
    virtual void setTargetValue(Item* target, V value) const = 0;
};

// hide the rest of this block of code from Doxygen
/// \cond

// Specialization for item list properties
template<> class TypedPropertyAccessor<vector<Item*>> : public PropertyAccessor
{
protected:
    TypedPropertyAccessor() {}

public:
    virtual vector<Item*> targetValue(const Item* target) const = 0;
    virtual void clearTargetValue(Item* target) const = 0;
    virtual void insertIntoTargetValue(Item* target, int index, Item* value) const = 0;
    virtual void removeFromTargetValue(Item* target, int index) const = 0;
};

/// \endcond
// end of hiding code from Doxygen

////////////////////////////////////////////////////////////////////

/** The TargetTypedPropertyAccessor template class sits at level two of the PropertyAccessor class
    hierarchy. It actually implements the property access interface for a value type V and for a
    target type T (which must inherit Item), both given as template parameters.
    TargetTypedPropertyAccessor inherits TypedPropertyAccessor with an appropriate value type. The
    template is specialized for several value types to allow special processing, such as casting
    pointers between Item and one of its subclasses, or translating an enumeration value to an
    integer and vice versa. Refer to the PropertyAccessor class for more information. */
template<typename V, class T, typename = void> class TargetTypedPropertyAccessor : public TypedPropertyAccessor<V>
{
    using Getter = V (T::*)() const;
    using Setter = void (T::*)(V);
    Getter _getter;
    Setter _setter;

public:
    TargetTypedPropertyAccessor(Getter getter, Setter setter)
    {
        _getter = getter;
        _setter = setter;
    }

    V targetValue(const Item* target) const override
    {
        auto tg = dynamic_cast<const T*>(target);
        return tg ? (tg->*_getter)() : V{};
    }
    void setTargetValue(Item* target, V value) const override
    {
        auto tg = dynamic_cast<T*>(target);
        if (tg) (tg->*_setter)(value);
    }
};

// hide the rest of this block of code from Doxygen
/// \cond

// Specialization for enumeration properties, so that we can translate enumeration values to integers, and vice versa
// The third template argument ensures that this specialization is selected only for enumeration types
template<typename V, class T>
class TargetTypedPropertyAccessor<V, T, std::enable_if_t<std::is_enum<V>::value>> : public TypedPropertyAccessor<int>
{
    using Getter = V (T::*)() const;
    using Setter = void (T::*)(V);
    Getter _getter;
    Setter _setter;

public:
    TargetTypedPropertyAccessor(Getter getter, Setter setter)
    {
        _getter = getter;
        _setter = setter;
    }

    int targetValue(const Item* target) const override
    {
        auto tg = dynamic_cast<const T*>(target);
        return tg ? static_cast<int>((tg->*_getter)()) : 0;
    }
    void setTargetValue(Item* target, int value) const override
    {
        auto tg = dynamic_cast<T*>(target);
        if (tg) (tg->*_setter)(static_cast<V>(value));
    }
};

// Specialization for value list types, which use a const reference in the getter return type
template<typename V, class T> class TargetTypedPropertyAccessor<vector<V>, T> : public TypedPropertyAccessor<vector<V>>
{
    using Getter = const vector<V>& (T::*)() const;
    using Setter = void (T::*)(vector<V>);
    Getter _getter;
    Setter _setter;

public:
    TargetTypedPropertyAccessor(Getter getter, Setter setter)
    {
        _getter = getter;
        _setter = setter;
    }

    vector<V> targetValue(const Item* target) const override
    {
        auto tg = dynamic_cast<const T*>(target);
        return tg ? (tg->*_getter)() : vector<V>{};
    }
    void setTargetValue(Item* target, vector<V> value) const override
    {
        auto tg = dynamic_cast<T*>(target);
        if (tg) (tg->*_setter)(value);
    }
};

// Specialization for item properties, so that we can cast the Item subclass instances to Item, and vice versa
template<typename V, class T> class TargetTypedPropertyAccessor<V*, T> : public TypedPropertyAccessor<Item*>
{
    using Getter = V* (T::*)() const;
    using Setter = void (T::*)(V*);
    Getter _getter;
    Setter _setter;

public:
    TargetTypedPropertyAccessor(Getter getter, Setter setter)
    {
        _getter = getter;
        _setter = setter;
    }

    Item* targetValue(const Item* target) const override
    {
        auto tg = dynamic_cast<const T*>(target);
        return tg ? (tg->*_getter)() : nullptr;
    }
    void setTargetValue(Item* target, Item* value) const override
    {
        auto val = dynamic_cast<V*>(value);
        auto tg = dynamic_cast<T*>(target);
        if (tg) (tg->*_setter)(val);
    }
};

// Specialization for item list properties, to implement the special setters,
// and so that we can cast the Item subclass instances to Item, and vice versa
template<typename V, class T>
class TargetTypedPropertyAccessor<vector<V*>, T> : public TypedPropertyAccessor<vector<Item*>>
{
    using Getter = const vector<V*>& (T::*)() const;
    using Clearer = void (T::*)();
    using Inserter = void (T::*)(int, V*);
    using Remover = void (T::*)(int);
    Getter _getter;
    Clearer _clearer;
    Inserter _inserter;
    Remover _remover;

public:
    TargetTypedPropertyAccessor(Getter getter, Clearer clearer, Inserter inserter, Remover remover)
    {
        _getter = getter;
        _clearer = clearer;
        _inserter = inserter;
        _remover = remover;
    }

    vector<Item*> targetValue(const Item* target) const override
    {
        vector<Item*> result;
        auto tg = dynamic_cast<const T*>(target);
        if (tg)
            for (const auto& it : (tg->*_getter)())
                if (it) result.push_back(it);
        return result;
    }
    virtual void clearTargetValue(Item* target) const override
    {
        auto tg = dynamic_cast<T*>(target);
        if (tg) (tg->*_clearer)();
    }
    virtual void insertIntoTargetValue(Item* target, int index, Item* value) const override
    {
        auto val = dynamic_cast<V*>(value);
        auto tg = dynamic_cast<T*>(target);
        if (val && tg) (tg->*_inserter)(index, val);
    }
    virtual void removeFromTargetValue(Item* target, int index) const override
    {
        auto tg = dynamic_cast<T*>(target);
        if (tg) (tg->*_remover)(index);
    }
};

/// \endcond
// end of hiding code from Doxygen

////////////////////////////////////////////////////////////////////

#endif

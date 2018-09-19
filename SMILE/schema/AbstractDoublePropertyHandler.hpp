/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ABSTRACTDOUBLEPROPERTYHANDLER_HPP
#define ABSTRACTDOUBLEPROPERTYHANDLER_HPP

#include "PropertyHandler.hpp"

////////////////////////////////////////////////////////////////////

/** AbstractDoublePropertyHandler serves as an abstract base class for handling SMILE data item
    properties of floating point types with an optional unit specification. It offers common
    functionality for derived classes supporting a single floating point value
    (DoublePropertyHandler) or a list of values (DoubleListPropertyHandler). */
class AbstractDoublePropertyHandler : public PropertyHandler
{
    // ================== Constructing & Destructing ==================

public:
    /** Constructs a property handler for the specified target item and property, with a given
        schema definition. */
    using PropertyHandler::PropertyHandler;

    // ================== Specific functions for this property type ==================

private:
    /** Returns the minimum value for the handled property. If no minimum value is specified in the
        schema definition, the function returns negative infinity. */
    double minValue() const;

    /** Returns the maximum value for the handled property. If no maximum value is specified in the
        schema definition, the function returns positive infinity. */
    double maxValue() const;

    /** Returns true if the minimum value for the handled property is excluded from the permitted
        range (indicated by a leading "]"), or if no minimum value is specified in the schema
        definition. */
    bool minOpen() const;

    /** Returns true if the maximum value for the handled property is excluded from the permitted
        range (indicated by a trailing "["), or if no maximum value is specified in the schema
        definition. */
    bool maxOpen() const;

public:
    /** Returns a human-readable description of the permitted value range for the handled property.
        The description includes brackets indicating an open or closed range on either side, and
        default unit specifications for physical quantities. If no minimum/maximum value is
        specified in the schema definition, the description uses negative/positive infinity. */
    string rangeDescription() const;

    /** Returns true if the given double value is finite and within the permitted value range for
        the handled property, and false otherwise. */
    bool isInRange(double value) const;

    /** Returns true if all values in the given double list are finite and within the permitted
        value range for the handled property, and false otherwise. */
    bool isInRange(const vector<double>& value) const;

    /** Returns true if the specified string is non-empty and contains a valid string
        representation of a floating point number with an optional unit specification. Otherwise
        returns false. If present, the unit specification must follow the number, separated by one
        or more spaces. The allowed set of unit specifications is derived from the physical
        quantity attribute of the handled property; the default unit specification is determined
        from the unit system associated with the simulation hierarchy in which the handled property
        resides. */
    bool isValidDouble(string value) const;

    /** Returns the double value represented by the specified string, or zero if the string is
        empty or contains an invalid representation. See isValid() for more information. */
    double toDouble(string value) const;

    /** Returns a string representation of the specified double value, including an appropriate
        unit specification. See isValid() for more information. */
    string toString(double value) const;

    /** Returns true if the specified string is non-empty and contains a comma separated list, with
        each item in the list representing a valid floating point number with an optional unit
        specification, according to the format described for the isValidDouble() function.
        Otherwise returns false. */
    bool isValidDoubleList(string value) const;

    /** Returns the list of double values represented by the specified string. If the string is
        empty, or any of the comma-separated items in the specified string are invalid, the
        function returns an empty list. See isValidDoubleList() for more information. */
    vector<double> toDoubleList(string value) const;

    /** Returns a string representation of the specified list of double values, each item including
        an appropriate unit specification, and commas seperating the items. See isValidDoubleList()
        for more information. */
    string toString(vector<double> value) const;

    /** Returns the physical quantity name for the handled property, as described in the schema
        definition class, or the empty string if the handled property is a dimensionless quantity.
        If the value of the "quantity" attribute in the schema definition starts with an at sign,
        the quantity is determined instead as the string value of the indicated enumeration
        property. In that case, because the string value of an enumeration property can't be empty,
        a value representing an unknown quantity is silently replaced by the empty string to
        represent a dimensionless quantity. */
    string quantity() const;

    /** Returns the name of the unit system associated with the dataset in which the handled
        property resides. If the schema definition does not provide a unit system, the function
        throws an error. If the schema definition provides a single unit system, the function
        returns the name of that unit system. If the schema definition provides two or more unit
        systems, the function searches the hierarchy of the target dataset for a SMILE data item
        that inherits the common unit system base type, and returns the actual type of that item.
        If no such data item can be found, the function throws an error. */
    string unitSystem() const;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef UNITDEF_HPP
#define UNITDEF_HPP

#include "Basics.hpp"
#include <map>

////////////////////////////////////////////////////////////////////

/** The UnitDef class represents information about the units and unit systems optionally used in a
    SMILE data set. It offers functions for adding unit and unit system information during its
    construction, and for querying the accumulated information in various ways.

    The information held by a UnitDef instance (commonly called a unit definition) can be
    summarized as follows. There is a list of physical quantities. Each physical quantity has one
    or more units associated with it. Values in these units can be converted to and from internal
    program units through a linear transformation. Furthermore, there is a list of one or more unit
    systems. A unit system defines the default unit for each physical quantity.

    The UnitDef class serves two distinct (but related) purposes. The SchemaDef class uses a
    UnitDef instance to store information about the units and unit systems used in the SMILE data
    set being described by the schema definition. The SchemaDef class has \em friend status so that
    it can freely access all unit definition contents when loading and saving schema definitions
    from and to XML files.

    In addition, the UnitDef class plays an important role when SMILE library client code wishes to
    define and support units as part of a hierarchy of client SMILE subclasses (refer to the
    documentation of the Item class for information on client items). In this case, the client code
    developer should proceed as follows.

    - Define a UnitDef subclass with a constructor that loads the unit and unit system information.
    - Register this subclass with the item registry using the ItemRegistry::addUnitDef() function.
    - Define an Item subclass with the same name as each of the unit systems loaded into the unit
      definition, with a common "unit system" base class (in between Item and these classes).
    - Include these unit system classes in the overall Item type hierarchy, allowing the presence
      of one of these types in a SMILE dataset to determine the unit system used in that dataset.

    If only a single unit system was loaded into the unit definition, the last two steps are
    not required (but still allowed), because that unique unit system will be used by default.

    \attention <i>All names used in this class for physical quantities, units, and unit systems
    are <b>case sensitive</b> and are <b>not allowed to contain white space</b>.</i>

    Here is an example of a UnitDef subclass constructor loading unit definition information.

    \verbatim
    addUnit("wavelength", "m", 1.);
    addUnit("wavelength", "micron", 1e-6);
    addUnit("wavelength", "Angstrom", 1e-10);
    addUnit("temperature", "K", 1.);
    addUnit("temperature", "C", 1., 1., 273.15);

    addDefaultUnit("SIUnits", "wavelength", "m");
    addDefaultUnit("SIUnits", "temperature", "K");
    addDefaultUnit("StellarUnits", "wavelength", "micron");
    addDefaultUnit("StellarUnits", "temperature", "K");
    \endverbatim
*/
class UnitDef
{
    friend class SchemaDef;  // so that schema definition can load and extract info at will

    // ================== Constructing and assigning ==================

public:
    /** This constructor creates an empty unit definition. Use the addUnit() and addDefaultUnit()
        functions from the constructor of a subclass to load the appropriate information into the
        definition. */
    UnitDef();

    /** The copy constructor creates a unit definition that contains a copy of the information in
        the given unit definition. */
    UnitDef(const UnitDef& UnitDef) = default;

    /** The assignment operator copies the information from the given unit definition, discarding
        any information that may have pre-existed in the receiving definition. */
    UnitDef& operator=(const UnitDef&) = default;

    // ================== Loading information ==================

protected:
    /** This function adds information about a particular unit. It can be called from the
        constructor of a subclass to load the appropriate information into the definition. In order
        of occurrence, the arguments specify the name of the physical quantity using this unit, the
        name of the unit, and the factor \f$f\f$, power index \f$p\f$ and offset \f$o\f$ needed to
        convert a value in this unit to a value in the corresponding internal program unit. The
        conversion is performed using \f[ v_\mathrm{program} = f\,v^p + o.\f] If not specified, the
        power index defaults to one and the offset defaults to zero. */
    void addUnit(string quantity, string unit, double factor, double power = 1., double offset = 0.);

    /** This function specifies the default unit for a particular quantity in a given unit system.
        It can be called from the constructor of a subclass to load the appropriate information
        into the definition. In order of occurrence, the arguments specify the name of the unit
        system, the name of the physical quantity, and the name of the corresponding default unit.
        */
    void addDefaultUnit(string unitSystem, string quantity, string unit);

    // ================== Retrieving information ==================

public:
    /** This function returns returns true if the specified physical quantity is present in the
        unit definition, and false if it is not. */
    bool has(string qty) const;

    /** This function returns true if the specified combination of physical quantity and unit or
        unit system is present in the unit definition, and false if not. The name of the physical
        quantity must always be specified. The unit can be specified either directly, or
        indirectly by providing the name of a unit system. In the latter case, the function uses
        the default unit for the specified quantity in the specified unit system. */
    bool has(string qty, string unit) const;

    /** This function returns the definition of the specified combination of physical quantity and
        unit or unit system in the form of a tuple providing the front factor, power exponent and
        offset for conversion from input to internal quantities. Refer to the has() function for a
        description of how to specify the physical quantity and unit. If the specified combination
        is not present in the unit definition, the function throws an exception. */
    std::tuple<double, double, double> def(string qty, string unit) const;

    /** This function converts a physical value from the specified units to internal program units.
        Refer to the has() function for a description of how to specify the physical quantity and
        unit. If the specified combination is not present in the unit definition, the function
        throws an exception. */
    double in(string qty, string unit, double value) const;

    /** This function converts a physical value from internal program units to the specified units.
        Refer to the has() function for a description of how to specify the physical quantity and
        unit. If the specified combination is not present in the unit definition, the function
        throws an exception. */
    double out(string qty, string unit, double value) const;

    /** This function returns the name of the default unit listed in the schema definition for the
        specified physical quantity in the specified unit system. If the specified combination of
        physical quantity and unit system is not present in the unit definition, the function
        throws an exception. */
    string unit(string qty, string unitSystem) const;

    // ================== Data members ==================

private:
    // use ordered maps below so that things are sorted when retrieved by the schema definition for saving
    // (the implementation of the UnitDef class does not depend on the maps being ordered)

    // a list of all physical quantities with corresponding units:
    //      <quantity-name, <unit-name, (factor, power, offset)>>
    std::map<string, std::map<string, std::tuple<double, double, double>>> _quantities;

    // a list of all units systems with corresponding default units:
    //      <unitsystem-name, <quantity-name, unit-name>>
    std::map<string, std::map<string, string>> _unitSystems;
};

////////////////////////////////////////////////////////////////////

#endif

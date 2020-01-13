/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SCHEMADEF_HPP
#define SCHEMADEF_HPP

#include "Basics.hpp"
#include "TypeDef.hpp"
#include "UnitDef.hpp"
#include <map>
class Item;
class NameManager;
class PropertyHandler;

////////////////////////////////////////////////////////////////////

/** The SchemaDef class represents a SMILE schema definition. It offers facilities for loading a
    SMILE schema from file or through calling the relevant "add" functions, and for querying the
    information in the loaded schema in various ways. */
class SchemaDef final
{
    // ================== Getting info on a schema ==================

public:
    /** This static function returns the value of the 'title' attribute on the Schema element in
        the SMILE XML file with the specified file path, or the empty string if an error occurs. */
    static string getSchemaTitle(string filePath);

    /** This static function determines whether the SMILE dataset with the specified file path is
        compatible with the SMILE schema defined in the SMILE XML file with the specified file
        path, by examining the root and top-level elements in both XML files. The function returns
        true if the dataset is compatible with the schema, and false if not or if an error occurs.
        */
    static bool isCompatible(string schemaFilePath, string dataFilePath);

    // ================== Loading and saving a schema ==================

public:
    /** This constructor loads a particular SMILE schema from the SMILE XML file with the specified
        path. When an error occurs while opening, parsing or interpreting the XML data stream, the
        constructor throws a FatalError with an appropriate error message. */
    SchemaDef(string filePath);

    /** This constructor creates an empty SMILE schema with the specified basic properties. The
        type, property and and unit information for the schema must be added by calling the
        appropriate "add" functions. The constructor arguments specify, in order of occurrence:
        - \em name: a short name for this schema
        - \em title: a description of the type of datasets described by this schema
        - \em version: the version of this schema definition
        - \em extension: the filename extension for datasets described by this schema
        - \em root: the name of the root element in datasets described by this schema
        - \em type: the type of the top-level item in datasets described by this schema
        - \em format: the version of the described data format (specified on the root element)
        - \em url: a URL pointing to information on the Web for this schema (or the empty string if not available)
    */
    SchemaDef(string name, string title, string version, string extension, string root, string type, string format,
              string url);

    /** This function creates an empty type definition, adds it to the schema definition, and
        returns a writable reference to it. The first three arguments specify the name of the type,
        the name of its base type, and the description of the type. The last argument specifies the
        function that will be used to create instances of this type. It should be provided for all
        concrete types, and should be left to the default (a null pointer) for abstract types. The
        type name, the base type name, the description and the concrete/abstract flag are set by
        this function. All other information in the type definition, including any property
        definitions, must be completed or added by calling the type definition's setter functions.
        This function throws a fatal error if the schema definition already contains a type with
        the specified name (to keep the function from being abused to retrieve a writable reference
        to an existing type definition). */
    TypeDef& addTypeDef(string name, string base, string title, TypeDef::Instantiator instantiator = nullptr);

    /** This function copies the contents of the specified unit definition into the schema
        definition. Any pre-existing unit and unit system information is discarded. */
    void loadUnitDef(const UnitDef& unitDef);

    /** This function saves the SMILE schema to the SMILE XML file with the specified path. The
        optional last argument specifies a producer identification string to be included as an
        attribute on the root element. When an error occurs while opening the file or while writing
        the XML data stream, the function throws a FatalError with an appropriate error message. */
    void save(string filePath, string producer = string()) const;

    // ================== Copying and moving ==================

public:
    /** The copy constructor is deleted. SchemaDef instances can't be copied or moved because of
        the PropertyDef instances they contain. */
    SchemaDef(const SchemaDef&) = delete;

    /** The assignment operator is deleted. SchemaDef instances can't be copied or moved because of
        the PropertyDef instances they contain. */
    SchemaDef& operator=(const SchemaDef&) = delete;

    // ================== Retrieving schema information ==================

public:
    /** Returns a description of the producer of the file from which this schema was loaded. */
    string schemaProducer() const;

    /** Returns a short name for this schema. */
    string schemaName() const;

    /** Returns a description of the type of datasets described by this schema. */
    string schemaTitle() const;

    /** Returns the version of this schema definition. */
    string schemaVersion() const;

    /** Returns the filename extension for datasets described by this schema. */
    string schemaExtension() const;

    /** Returns the name of the root element in datasets described by this schema. */
    string schemaRoot() const;

    /** Returns the type of the top-level item in datasets described by this schema. */
    string schemaType() const;

    /** Returns the version of the described data format (specified on the root element in data
        sets described by this schema). */
    string schemaFormat() const;

    /** Returns a URL pointing to information on the Web for this schema. */
    string schemaUrl() const;

    // -------------------------------------

    /** Returns the title (used for display to a user) associated with the specified type. The
        function throws an error if the specified type is not defined in the schema. */
    string title(string type) const;

    /** Returns the titles (used for display to a user) associated with the specified types, in the
        same order. */
    vector<string> titles(const vector<string>& types) const;

    /** Returns true if the first type inherits the second. The function throws an error if the
        first type (the child type) is not defined in the schema. */
    bool inherits(string childType, string parentType) const;

    /** Returns a list of the types from which the specified type inherits, directly or indirectly,
        starting with the type itself up to and including the root type of the hierarchy. The list
        includes both concrete and abstract types. The function throws an error if the specified
        type is not defined in the schema. */
    vector<string> ascendants(string type) const;

    /** Returns a list of concrete types that inherit the specified type, in the order listed in
        the schema definition. The function throws an error if the specified type is not defined in
        the schema. */
    vector<string> descendants(string type) const;

    /** Returns the names of all properties for the specified type, including inherited properties
        for all direct and indirect base types. By default, base type properties are listed first,
        and for each type, properties are listed in the order of the schema definition. However,
        base types can specify that subtype properties must be inserted before a given property
        rather than being added at the end of the list. */
    vector<string> properties(string type) const;

    /** Returns the name of the particular (base) type in which the specified property is defined
        for the specified (concrete) type. If the specified property is not defined in the
        specified type or in one of its base types, the function throws a fatal error. */
    string definingType(string type, string property) const;

    /** Returns the title (used for display to a user) associated with the specified property of
        the specified type. The function throws an error if the specified property and type
        combination is not defined in the schema. */
    string propertyTitle(string type, string property) const;

    // -------------------------------------

    /** Returns a Boolean expression that, when evaluated against the current global and local name
        sets, will determine whether the specified type is allowed. To obtain this result, the
        Boolean expressions in the "allowedIf" attribute values for the specified type and for any
        of its base types, recursively, are concatenated with the AND operator into a single
        Boolean expression. The function throws an error if the specified type is not defined in
        the schema. */
    string allowed(string type) const;

    /** Returns a Boolean expression that, when evaluated against the current global and local name
        sets, will determine whether the specified type is allowed and displayed. To obtain this
        result, the Boolean expressions in the "allowedIf" and "displayedIf" attribute values for
        the specified type and for any of its base types, recursively, are concatenated with the
        AND operator into a single Boolean expression. The function throws an error if the
        specified type is not defined in the schema. */
    string allowedAndDisplayed(string type) const;

    /** Returns a list of conditional value expressions that, when evaluated against the current
        global and local name sets, will provide the names that need to be inserted when the
        specified type is entered into a dataset. The list includes the conditional value
        expressions given as "insert" attribute values for the specified type and for any of its
        base types, recursively. Empty "insert" attribute values are not included, so the returned
        list may be empty. The function throws an error if the specified type is not defined in the
        schema. */
    vector<string> toBeInserted(string type) const;

    // -------------------------------------

    /** Creates a new SMILE data item of the specified type. Ownership of the new object is passed
        to the caller through a unique pointer. This function throws an error if the specified type
        is not defined in the schema, or if the type is an abstract type. */
    std::unique_ptr<Item> createItem(string type) const;

    /** Creates a new property handler for the specified property in the specified SMILE data item,
        and returns a unique pointer to this new instance. The third argument, a pointer to the
        name manager used for the current dataset, is passed to the property handler so that it can
        access the global and local name sets for evaluating Boolean expressions.

        The returned handler is of a PropertyHandler subclass appropriate for the property type.
        Because it is guarded by a unique pointer, the handler is automatically deleted when the
        return value goes out of scope. This function throws an error if the specified item does
        not have a property with the specified name. */
    std::unique_ptr<PropertyHandler> createPropertyHandler(Item* item, string property, NameManager* nameMgr) const;

    // -------------------------------------

    /** This function returns returns true if the specified physical quantity is provided in the
        schema definition, and false if it is not. The name of the physical quantity is case
        sensitive and should not contain any spaces. */
    bool has(string qty) const;

    /** This function returns true if the specified combination of physical quantity and unit or
        unit system is provided in the schema definition, and false if not. The name of the
        physical quantity must always be specified. The unit can be specified either directly, or
        indirectly by providing the name of a unit system. In the latter case, the function uses
        the default units for the specified quantity in the specified unit system. All names
        (physical quantity, unit system, and units) are case sensitive and should not contain any
        spaces. */
    bool has(string qty, string unit) const;

    /** This function converts a physical value from the specified units to internal program units.
        Refer to the has() function for a description of how to specify the physical quantity and
        unit. If the specified combination is not provided in the schema definition, the function
        throws an exception. */
    double in(string qty, string unit, double value) const;

    /** This function converts a physical value from internal program units to the specified units.
        Refer to the has() function for a description of how to specify the physical quantity and
        unit. If the specified combination is not provided in the schema definition, the function
        throws an exception. */
    double out(string qty, string unit, double value) const;

    /** This function returns the name of the default unit listed in the schema definition for the
        specified physical quantity in the specified unit system. The physical quantity and unit
        system names are case sensitive and should not contain any spaces. If the specified
        combination of physical quantity and unit system is not provided in the schema definition,
        the function throws an exception. */
    string unit(string qty, string unitSystem) const;

    /** This function returns true if two or more unit systems are provided in the schema
        definition. Otherwise it returns false. */
    bool hasMultipleUnitSystems() const;

    /** If the schema definition does not provide a unit system, this function throws an error. If
        the schema definition provides a single unit system, the function returns the name of that
        unit system. If the schema definition provides two or more unit systems, the function
        returns the name of the common base type for the types corresponding to each of the unit
        systems. If not all of the unit systems have corresponding types, or if there is no common
        base type, the function throws an error. The common base type is defined as a type with a
        set of concrete descendents that exactly matches the set of unit systems. */
    string unitSystemBase() const;

    // ================== Private utilities ==================

private:
    /** Returns a reference to the type definition for the specified type, or throws an
        error if the specified type is not defined in the schema. */
    const TypeDef& typeDef(string type) const;

    /** Returns a reference to the property definition for the specified type and property, or
        throws an error if the specified type and property combination is not defined in the
        schema. */
    const PropertyDef& propertyDef(string type, string property) const;

    // ================== Data members ==================

private:
    string _producer;   // a description of the producer of the file from which this schema was loaded
    string _name;       // a short name for this schema
    string _title;      // a description of the type of datasets described by this schema
    string _version;    // the version of this schema definition
    string _extension;  // the filename extension for datasets described by this schema
    string _root;       // the name of the root element in datasets described by this schema
    string _type;       // the type of the top-level item in datasets described by this schema
    string _format;     // the version of the described data format (specified on the root element)
    string _url;        // a URL pointing to information on the Web for this schema (or the empty string)

    // use ordered map so that things are sorted in XML files generated by save() function
    std::map<string, TypeDef> _allTypes;  // a list of all types including concrete and abstract types
    vector<string> _concreteTypes;        // a list of all concrete types in order of definition

    UnitDef _unitDef;  // the information on units and unit systems held by this schema definition
};

////////////////////////////////////////////////////////////////////

#endif

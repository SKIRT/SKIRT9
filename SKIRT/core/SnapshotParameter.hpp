/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SNAPSHOTPARAMETER_HPP
#define SNAPSHOTPARAMETER_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** SnapshotParameter is a helper class for holding metadata about a configurable parameter to be
    imported from a snapshot text column file. For example, an SEDFamily object would return a list
    of SnapshotParameter objects to define the parameters needed to select a particular %SED
    template from the family. Similarly, a MaterialMix object can request parameters to be imported
    from a snapshot.

    The SnapshotParameter class offers specific support for a number of frequently-used snapshot
    parameters, so that these can be discovered and, for example, used for probing. Examples
    include metallicity and age. Use one of the factory functions to construct a SnapshotParameter
    instance. */
class SnapshotParameter
{
    // ================== Data types ==================

public:
    /** This enumeration lists the identifiers for the supported snapshot parameter types. */
    enum class Identifier { InitialMass, CurrentMass, Metallicity, Age, Temperature, Custom };

    // ================== Constructing ==================

private:
    /** This constructor initializes a snapshot parameter object with the given information for
        each field. The constructor is private; use one of the factory functions instead. */
    SnapshotParameter(Identifier identifier, string description, string quantity, string defaultUnit);

public:
    /** This function returns a SnapshotParameter instance of type InitialMass with default units
        of solar mass. */
    static SnapshotParameter initialMass();

    /** This function returns a SnapshotParameter instance of type CurrentMass with default units
        of solar mass. */
    static SnapshotParameter currentMass();

    /** This function returns a SnapshotParameter instance of type Metallicity. */
    static SnapshotParameter metallicity();

    /** This function returns a SnapshotParameter instance of type Age with default units of year.
        */
    static SnapshotParameter age();

    /** This function returns a SnapshotParameter instance of type Temperature with default units
        of Kelvin. */
    static SnapshotParameter temperature();

    /** This function returns a custom SnapshotParameter instance. The \em description argument is
        used only for logging purposes. The \em quantity argument specifies the physical quantity
        represented by the parameter. It must match one of the quantity strings supported by the
        Units system, or one of the special quantity strings recognized by the
        TextInFile::addColumn() function. The \em defaultUnit argument specifies the default unit
        string, which is used in case the input file does not contain unit information for the
        parameter. */
    static SnapshotParameter custom(string description, string quantity = string(), string defaultUnit = string());

    // ================== Querying ==================

public:
    /** This function returns the identifier for the snapshot parameter. All custom snapshot
        parameters have the same identifier; there is no formal way to tell them apart. */
    Identifier identifier() const { return _identifier; }

    /** This function returns the description of the snapshot parameter. */
    const string& description() const { return _description; }

    /** This function returns the string defining the physical quantity represented by the snapshot
        parameter. */
    const string& quantity() const { return _quantity; }

    /** This function returns the default units to be used for the snapshot parameter in case the
        input file does not contain unit information for it. */
    const string& defaultUnit() const { return _defaultUnit; }

    // ================== Data Members ==================

private:
    Identifier _identifier;
    string _description;
    string _quantity;
    string _defaultUnit;
};

////////////////////////////////////////////////////////////////////

#endif

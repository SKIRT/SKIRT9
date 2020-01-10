/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SNAPSHOTPARAMETER_HPP
#define SNAPSHOTPARAMETER_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** SnapshotParameter is an almost trivial helper class for holding metadata about a configurable
    parameter to be imported from a snapshot text column file. For example, an SEDFamily object
    would return a list of SnapshotParameter objects to define the parameters needed to select a
    particular SED template from the family. */
class SnapshotParameter
{
public:
    /** This constructor initializes a snapshot parameter object with the given description and
        unit information. The \em description argument is used only for logging purposes. The \em
        quantity argument specifies the physical quantity represented by the parameter. It must
        match one of the quantity strings supported by the Units system, or one of the special
        quantity strings recognized by the TextInFile::addColumn() function. The \em defaultUnit
        argument specifies the default unit string, which is used in case the input file does not
        contain unit information for the parameter. */
    SnapshotParameter(string description, string quantity = string(), string defaultUnit = string())
        : _description(description), _quantity(quantity), _defaultUnit(defaultUnit)
    {}

    /** This function returns the parameter description. */
    const string& description() const { return _description; }

    /** This function returns the string defining the physical quantity represented by the
        parameter. */
    const string& quantity() const { return _quantity; }

    /** This function returns the default units to be used for the parameter in case the input file
        does not contain unit information for the parameter. */
    const string& defaultUnit() const { return _defaultUnit; }

private:
    string _description;
    string _quantity;
    string _defaultUnit;
};

////////////////////////////////////////////////////////////////////

#endif

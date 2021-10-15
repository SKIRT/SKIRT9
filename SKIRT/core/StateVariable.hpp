/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATEVARIABLE_HPP
#define STATEVARIABLE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** StateVariable is a helper class for holding metadata about a medium state variable. For
    example, the MaterialMix::specificStateVariableInfo() function returns a list of StateVariable
    instances describing the specific state variables used by the receiving material mix. This
    allows the MediumSystem class to allocate storage for the appropriate set of state variables,
    and it allows probing the relevant medium state variables for output.

    For each spatial cell in the simulation, the medium state includes a set of \em common state
    variables shared by all medium components and a set of \em specific state variables for each
    individual medium component. Refer to the MediumState class for a list of state variables that
    may be present in each of these sets.

    Each StateVariable instance includes an identifier specifying one of the supported state
    variables. Multiple custom variables can be specified by adding an index \f$0 \le k < K\f$ to
    the base identifier, where \f$K\f$ is the total number of custom variables. The StateVariable
    instance also includes a human-readable description of the quantity being represented, a unit
    system quantity definition (string), and a text output format specifier. This information
    allows custom variables to be output by generic probes or to be handled by portions of the code
    that cooperate specifically with this type of material mix, such as recipes for adjusting the
    medium state based on the calculated radiation field.

    This class includes explicit support for the known state variables. Use one of the factory
    functions to construct a StateVariable instance. */
class StateVariable
{
    // ================== Data types ==================

public:
    /** This enumeration lists the identifiers for the supported state variables as indicated in
        the table in the class header. */
    enum class Identifier { Volume, BulkVelocity, MagneticField, NumberDensity, Metallicity, Temperature, Custom };

    // ================== Constructing ==================

private:
    /** This constructor initializes a state variiable object with the given information for each
        field. The constructor is private; to construct a StateVariable instance, use one of
        the factory functions instead. */
    StateVariable(Identifier identifier, int customIndex, string description, string quantity, char format);

public:
    /** This function returns a StateVariable instance of type Volume. */
    static StateVariable volume();

    /** This function returns a StateVariable instance of type BulkVelocity. */
    static StateVariable bulkVelocity();

    /** This function returns a StateVariable instance of type MagneticField. */
    static StateVariable magneticField();

    /** This function returns a StateVariable instance of type NumberDensity. */
    static StateVariable numberDensity();

    /** This function returns a StateVariable instance of type Metallicity. */
    static StateVariable metallicity();

    /** This function returns a StateVariable instance of type Temperature. */
    static StateVariable temperature();

    /** This function returns a StateVariable instance of type Custom with the specified index and
        descriptive information. The latter includes a human-readable description of the quantity
        being represented, a unit system quantity definition, and a text output format specifier,
        i.e. 'd' for integer notation (even if the value is stored as a double), 'f' for fixed
        point notation, 'e' for scientific notation, and 'g' for the most concise 'f' or 'e'. The
        default format specifier is 'e'. */
    static StateVariable custom(int customIndex, string description, string quantity, char format = 'e');

    // ================== Querying ==================

public:
    /** This function returns the identifier for the state variable. All custom variables have the
        same identifier; they can be told apart through the customIndex() function(). */
    Identifier identifier() const { return _identifier; }

    /** This function returns the custom index for the state variable if it is of type Custom, or
        zero otherwise. */
    int customIndex() const { return _customIndex; }

    /** This function returns human-readable description for the state variable. */
    const string& description() const { return _description; }

    /** This function returns the string defining the physical quantity represented by the state
        variable. */
    const string& quantity() const { return _quantity; }

    /** This function returns the format specifier for the state variable, i.e. 'd' for integer
        notation (even if the value is stored as a double), 'f' for fixed point notation, 'e' for
        scientific notation, and 'g' for the most concise 'f' or 'e'. */
    char format() const { return _format; }

private:
    Identifier _identifier;
    int _customIndex;
    string _description;
    string _quantity;
    char _format;
};

////////////////////////////////////////////////////////////////////

#endif

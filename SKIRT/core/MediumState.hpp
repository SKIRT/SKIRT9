/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSTATE_HPP
#define MEDIUMSTATE_HPP

#include "Array.hpp"
#include "StateVariable.hpp"
#include "Vec.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the MediumState class contains the complete medium state for a simulation. The
    MediumSystem maintains a single (private) instance of this class and provides access to its
    contents to the rest of the simulation code as needed.

    For each spatial cell in the simulation, the medium state includes a set of \em common state
    variables shared by all medium components and a set of \em specific state variables for each
    individual medium component. The following table lists the variables that may be present in
    each of these sets.

    | Identifier | Symbol | Type | Set | Present
    | -----------|--------|------|-----|--------
    | Volume        | \f$V\f$      | double | Common   | always
    | BulkVelocity  | \f$v\f$      | Vec    | Common   | if needed for any medium component
    | MagneticField | \f$\bf{B}\f$ | Vec    | Common   | if needed for any medium component
    | NumberDensity | \f$n\f$      | double | Specific | always
    | Temperature   | \f$T\f$      | double | Specific | if requested by the medium component
    | Custom + k    |              | double | Specific | zero or more as requested by the medium component

    The presence or absence of these state variables is communicated to the MediumState class
    during its construction, so that it can dynamically allocate the appropriate amount of storage
    and provide an access mechanism for each state variable. This process is described in more
    detail in the following sections.

    <b>Construction</b>

    Fully initializing a MediumState object requires calling the following functions in the
    correct order (after default-constructing the object):
     - the initConfiguration() function specifies the number of spatial cells and the number of
       medium components in the simulation.
     - the initCommonStateVariables() function specifies the set of common state variables.
     - the initSpecificStateVariables() function must be called once for each medium component, in
       order of component index, specifying the set of specific state variables for that component.
     - the initAllocate() function finalizes construction and actually allocates storage.
     - the initCommunicate() function communicates the state variable values between processes.

    The initCommonStateVariables() and initSpecificStateVariables() functions each receive a list
    of StateVariable objects to identify the required state variables (only the identifier is used,
    the other information is ignored). All required state variables must be listed, including those
    that should always present. Variables of type Custom must be listed last. Multiple variables of
    type Custom can be requested by supplying indices in the range \f$ 0 \le k < K\f$, where K is
    the total number of custom variables. Each of these indices must occur in the list exactly once
    in increasing order.

    <b>Storage</b>

    To simplify the storage mechanism, state variables of type \c Vec are split into their three
    components, so that all state variables can considered to be of type \c double. Memory is
    allocated only for state variables that are actually needed. Given the number of spatial cells
    \f$M\f$, the number of medium components \f$H\f$, the number of common state variables \f$C\f$,
    and the number of specific state variables \f$S_h\f$ for each medium component \f$h\f$, the
    number of state variables per spatial cell is \f$K = C+\sum_h S_h\f$ and the grand total number
    of state variables is \f$N = M (C+\sum_h S_h)\f$.

    The MediumState class allocates a single one-dimensional data array of this size \f$N\f$ and
    provides a mapping to locate a particular state variable in this array. This is accomplished by
    also storing the offset for each variable within the set of variables per spatial cell. To
    facilitate access, these offsets are stored for all supported variables (by name) and not just
    for the required variables. The offsets for unused variables remain at zero, causing undefined
    behavior if client code attempts to access a variable that has not been requested.

    Given the offset \f$O_x\f$ for a particular state variable \f$x\f$, a spatial cell index
    \f$m\f$ and a medium component index \f$h\f$, the index of the variable in the data array can
    be calculated as \f$i=K \times m + O_x\f$ (storing contiguously per cell) or \f$i=M \times O_x
    + m\f$ (storing contiguously per variable). We can evaluate the performance of these and
    possibly other mapping schemes. */
class MediumState
{
    //============= Construction =============

public:
    /** This function initializes the number of spatial cells and number of medium components. */
    void initConfiguration(int numCells, int numComps);

    /** This function initializes the set of required common state variables. */
    void initCommonStateVariables(const vector<StateVariable>& variables);

    /** This function initializes the set of required specific state variables for the next medium
        component. The function must be called once for each medium component, in order of
        component index. */
    void initSpecificStateVariables(const vector<StateVariable>& variables);

    /** This function ends the initialization sequence, allocates memory for the state variables,
        and returns the total number \f$N = M (C+\sum_h S_h)\f$ of state variables. All newly
        allocated state variables are guaranteed to have a value of zero. */
    size_t initAllocate();

    /** This function communicates the state variable values between processes after each process
        has initialized the values for a subset of the spatial cells and left the values for the
        other cells at zero. (The function uses the ProcessManager::sumToAll() function, so it
        assumes that the uninitialized variables have a zero value). */
    void initCommunicate();

    //============= Setting =============

public:
    /** This function sets the volume \f$V\f$ of the spatial cell with index \f$m\f$. */
    void setVolume(double value, int m);

    /** This function sets the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in the
        spatial cell with index \f$m\f$. */
    void setBulkVelocity(Vec value, int m);

    /** This function sets the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell with index
        \f$m\f$. */
    void setMagneticField(Vec value, int m);

    /** This function sets the number density of the medium component with index \f$h\f$ in the
        spatial cell with index \f$m\f$. */
    void setNumberDensity(double value, int m, int h);

    /** This function sets the temperature \f$T\f$ of the medium component with index \f$h\f$ in
        the spatial cell with index \f$m\f$. */
    void setTemperature(double value, int m, int h);

    /** This function sets the value of the custom variable with index \f$i\f$ of the medium
        component with index \f$h\f$ in the spatial cell with index \f$m\f$. */
    void setCustom(double value, int i, int m, int h);

    //============= Querying =============

public:
    /** This function returns the volume \f$V\f$ of the spatial cell with index \f$m\f$. */
    double volume(int m) const { return _data[_numVars * m + _off_volu]; }

    /** This function returns the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in
        the spatial cell with index \f$m\f$. */
    Vec bulkVelocity(int m) const
    {
        int i = _numVars * m + _off_velo;
        return Vec(_data[i], _data[i + 1], _data[i + 2]);
    }

    /** This function returns the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell with
        index \f$m\f$. */
    Vec magneticField(int m) const
    {
        int i = _numVars * m + _off_mfld;
        return Vec(_data[i], _data[i + 1], _data[i + 2]);
    }

    /** This function returns the number density of the medium component with index \f$h\f$ in the
        spatial cell with index \f$m\f$. */
    double numberDensity(int m, int h) const { return _data[_numVars * m + _off_dens[h]]; }

    /** This function returns the temperature \f$T\f$ of the medium component with index \f$h\f$ in
        the spatial cell with index \f$m\f$. */
    double temperature(int m, int h) const { return _data[_numVars * m + _off_temp[h]]; }

    /** This function returns the value of the custom variable with index \f$i\f$ of the medium
        component with index \f$h\f$ in the spatial cell with index \f$m\f$. */
    double custom(int i, int m, int h) const { return _data[_numVars * m + _off_cust[h] + i]; }

    //======================== Data Members ========================

private:
    // data array containing the medium state variables
    Array _data;

    // configuration and offsets used for mapping to indices in the data array
    int _numCells{0};
    int _numComps{0};
    int _numVars{0};
    int _off_volu{0};
    int _off_velo{0};
    int _off_mfld{0};
    vector<int> _off_dens;
    vector<int> _off_temp;
    vector<int> _off_cust;

    // indices indicating the next item to be initialized; used only during initialization
    int _nextOffset{0};
    int _nextComponent{0};
};

////////////////////////////////////////////////////////////////

#endif

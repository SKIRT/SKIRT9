/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSTATE_HPP
#define MEDIUMSTATE_HPP

#include "Array.hpp"
#include "StateVariable.hpp"
#include "UpdateStatus.hpp"
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
    | Metallicity   | \f$Z\f$      | double | Specific | if requested by the medium component
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
     - the initAllocate() function finalizes construction and actually allocates storage;
       it initializes all variables to a value of zero.
     - the setXXX() functions set any nonzero initial variable values required to reflect the input
       model; this may happen in parallel.
     - the initCommunicate() function communicates the initialized state variable values between
       processes.

    The initCommonStateVariables() and initSpecificStateVariables() functions each receive a list
    of StateVariable objects to identify the required state variables. All required state variables
    must be listed, including those that should always present. Variables of type Custom must be
    listed last. Multiple variables of type Custom can be requested by supplying indices in the
    range \f$ 0 \le k < K\f$, where K is the total number of custom variables. Each of these
    indices must occur in the list exactly once in increasing order.

    <b>Synchronization</b>

    A simulation may update the values of medium state variables, often as part of an iterative
    process to calculate a self-consistent state. These updates may happen in parallel over cells.
    After a set of updates, usually at the end of each step in the iterative process, the
    synchronize() function should be called to broadcast the updates between processes.

    <b>Aggregation</b>

    When a simulation employs an iterative process to calculate a self-consistent medium state, it
    needs a criterion to determine whether the state has sufficiently converged. Such criterion is
    often based on aggregate information calculated over all cells in the spatial grid for the
    current iteration and for one or more previous iterations. This class can calculate and store
    such aggregate information in the form of "fake" cells that are added at the end of the array
    of regular cells. The number of stored aggregated states, including current and previous ones,
    must be specified as the last argument of the initConfiguration() function. If this number is
    zero, no aggregation is performed.

    Information from the aggregate states can be retrieved using the same functions as those for
    regular cells. Assuming a number of regular spatial cells \f$M\f$, the current (i.e. most
    recent) aggregate state has cell index \f$m=M\f$, the previous aggregate state has cell index
    \f$m=M+1\f$, the one before that has cell index \f$m=M+2\f$, and so forth.

    The calculateAggregate() function should be called at the end of construction, i.e. after
    initCommunicate(), and at the end of each dynamic medium state update cycle, i.e. after
    synchronize(). It recalculates the current aggregate state by accumulating information over all
    spatial cells. The pushAggregate() function should be called at the end of each step in the
    iterative process. It shifts the existing aggregate states to the next higher cell index,
    dropping the least recent aggregate state and making room for a new one.

    The current implementation performs the following aggregation. The cell volume is summed over
    all cells, providing the total volume of the spatial domain: \f$V_\mathrm{tot} = \sum V_m\f$.
    All variables with quantity type "numbervolumedensity" are summed after multiplication with the
    cell volume, providing a total number: \f$N_\mathrm{tot} = \sum n_m V_m\f$. Any other variables
    are not aggregated; i.e. the aggregated values are zero. The aggregated values can be retrieved
    through the volume(), numberDensity() and custom() functions using the fake cell index
    corresponding to the desired current or previous state.

    <b>Storage</b>

    To simplify the storage mechanism, state variables of type \c Vec are split into their three
    components, so that all state variables can considered to be of type \c double. Memory is
    allocated only for state variables that are actually needed. Given the number of spatial cells
    \f$M\f$ (possibly increased by the requested number of aggregate cells), the number of medium
    components \f$H\f$, the number of common state variables \f$C\f$, and the number of specific
    state variables \f$S_h\f$ for each medium component \f$h\f$, the number of state variables per
    spatial cell is \f$K = C+\sum_h S_h\f$ and the grand total number of state variables is \f$N =
    M (C+\sum_h S_h)\f$.

    The MediumState class allocates a single one-dimensional data array of this size \f$N\f$ and
    provides a mapping to locate a particular state variable in this array. This is accomplished by
    also storing the offset for each variable within the set of variables per spatial cell. To
    facilitate access, these offsets are stored for all supported variables (by name) and not just
    for the required variables. The offsets for unused variables remain at zero.

    Given the offset \f$O_x\f$ for a particular state variable \f$x\f$, a spatial cell index
    \f$m\f$ and a medium component index \f$h\f$, the index of the variable in the data array can
    be calculated as \f$i=K \times m + O_x\f$. This implies that variables are stored contiguously
    per cell.

    <b>Access to undefined variables</b>

    In general, an attempt to access (read or write) a variable for which storage has not been
    requested results in undefined behavior. However, assuming that storage has been requested for
    the volume (which is always needed anyway), the other variables in the common state can be
    safely read (not written) even if no storage was requested for them. This allows the bulk
    velocity and magnetic field vectors to be retrieved without concern for whether they were
    configured in the input model. */
class MediumState
{
    //============= Construction =============

public:
    /** This function initializes the number of spatial cells and number of medium components. If
        the specified number of aggregate cells is nonzero, the configuration is also prepared to
        store that number of aggregate states, as described in the class header. */
    void initConfiguration(int numCells, int numMedia, int numAggregateCells);

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

    //============= Synchronization and Aggregation =============

public:
    /** This function synchronizes the state variables between processes after each process has
        possibly updated some of their values. The provided vector indicates the cells for which
        the calling process has made some changes. The vector must be of length \em numCells as
        passed to the initConfiguration() function, where each element in the vector corresponds to
        the cell with the same index. If the calling process updated the state of a given cell, the
        corresponding vector element must be \em UpdatedConverged or \em UpdatedNotConverged and
        otherwise it must be \em NotUpdated. The function returns a pair of integers specifying the
        total number of updated cells and the total number of not-converged cells, aggregated over
        all processes.

        Only one of the calling processes may have updated the state for any given cell. If two or
        more processes updated the state of the same cell, the result of the synchronization is
        undefined. */
    std::pair<int, int> synchronize(const vector<UpdateStatus>& cellFlags);

    /** If aggregation has been requested when initializing the configuration, this function
        calculates the current aggregate state. For more information, see the description of
        aggregation in the class header. */
    void calculateAggregate();

    /** If aggregation has been requested when initializing the configuration, this function shifts
        the previously stored aggregate states to make room for a new current aggregate state. For
        more information, see the description of aggregation in the class header. */
    void pushAggregate();

    //============= Setting =============

public:
    /** This function sets the volume \f$V\f$ of the spatial cell with index \f$m\f$. */
    void setVolume(int m, double value);

    /** This function sets the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in the
        spatial cell with index \f$m\f$. */
    void setBulkVelocity(int m, Vec value);

    /** This function sets the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell with index
        \f$m\f$. */
    void setMagneticField(int m, Vec value);

    /** This function sets the number density of the medium component with index \f$h\f$ in the
        spatial cell with index \f$m\f$. */
    void setNumberDensity(int m, int h, double value);

    /** This function sets the metallicity \f$Z\f$ of the medium component with index \f$h\f$ in
        the spatial cell with index \f$m\f$. */
    void setMetallicity(int m, int h, double value);

    /** This function sets the temperature \f$T\f$ of the medium component with index \f$h\f$ in
        the spatial cell with index \f$m\f$. */
    void setTemperature(int m, int h, double value);

    /** This function sets the value of the custom variable with index \f$i\f$ of the medium
        component with index \f$h\f$ in the spatial cell with index \f$m\f$. */
    void setCustom(int m, int h, int i, double value);

    //============= Querying =============

public:
    /** This function returns the volume \f$V\f$ of the spatial cell with index \f$m\f$. */
    double volume(int m) const { return _data[_numVars * m + _off_volu]; }

    /** This function returns the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in
        the spatial cell with index \f$m\f$, or zero if storage was not requested for this
        variable. */
    Vec bulkVelocity(int m) const
    {
        if (_off_velo)
        {
            int i = _numVars * m + _off_velo;
            return Vec(_data[i], _data[i + 1], _data[i + 2]);
        }
        return Vec();
    }

    /** This function returns the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell with
        index \f$m\f$, or zero if storage was not requested for this variable. */
    Vec magneticField(int m) const
    {
        if (_off_mfld)
        {
            int i = _numVars * m + _off_mfld;
            return Vec(_data[i], _data[i + 1], _data[i + 2]);
        }
        return Vec();
    }

    /** This function returns the number density of the medium component with index \f$h\f$ in the
        spatial cell with index \f$m\f$. */
    double numberDensity(int m, int h) const { return _data[_numVars * m + _off_dens[h]]; }

    /** This function returns the metallicity \f$Z\f$ of the medium component with index \f$h\f$ in
        the spatial cell with index \f$m\f$. */
    double metallicity(int m, int h) const { return _data[_numVars * m + _off_meta[h]]; }

    /** This function returns the temperature \f$T\f$ of the medium component with index \f$h\f$ in
        the spatial cell with index \f$m\f$. */
    double temperature(int m, int h) const { return _data[_numVars * m + _off_temp[h]]; }

    /** This function returns the value of the custom variable with index \f$i\f$ of the medium
        component with index \f$h\f$ in the spatial cell with index \f$m\f$. */
    double custom(int m, int h, int i) const { return _data[_numVars * m + _off_cust[h] + i]; }

    //======================== Data Members ========================

private:
    // data array containing the medium state variables
    Array _data;

    // overall configuration
    int _numCells{0};
    int _numMedia{0};
    int _numAggregateCells{0};
    int _numVars{0};

    // offsets used for mapping common and specific variables (for each medium component) to indices in the data array
    int _off_volu{0};
    int _off_velo{0};
    int _off_mfld{0};
    vector<int> _off_dens;
    vector<int> _off_meta;
    vector<int> _off_temp;
    vector<int> _off_cust;

    // offsets used to aggregate all standard and custom specific variables of quantity type "numbervolumedensity"
    vector<int> _densityOffsets;

    // indices indicating the next item to be initialized; used only during initialization
    int _nextOffset{0};
    int _nextComponent{0};
};

////////////////////////////////////////////////////////////////

#endif

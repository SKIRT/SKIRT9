/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBEFORMBRIDGE_HPP
#define PROBEFORMBRIDGE_HPP

#include "Array.hpp"
#include "Direction.hpp"
#include "Position.hpp"
#include <functional>
class Form;
class Probe;
class SimulationItem;
class Snapshot;
class SpatialGrid;
class TextOutFile;
class Units;

//////////////////////////////////////////////////////////////////////

/** ProbeFormBridge is a helper class that provides a bridge between Probe and Form instances for
    form-assisted probes. The overall design objectives are described below. Refer to the
    documentation of the respective classes for more details.

    <b>Probes and forms</b>

    The idea is to allow combining two orthogonal concepts in an arbitrary fashion. The quantity
    being probed is specified by a Probe subclass, while the manner in which this quantity is being
    probed is specified by a Form subclass. Each Probe instance of this type is configured with a
    single Form instance.

    This design does not apply to all probes. It is intended for probing physical quantities that
    are defined as a scalar, vector or compound field over the spatial domain of the simulation.
    This includes two categories:

    - Spatial grid probes (inheriting the SpatialGridFormProbe class) handle quantities discretized
    over the simulation's spatial grid, including the radiation field and information stored in the
    medium state such as density, temperature, metallicity, velocity and custom variables. Because
    the quantities in this category are stored per spatial cell, they can all be probed by the same
    generic method.

    - Input model probes (inheriting the InputModelFormProbe class) handle quantities defined by
    the input model, provided the relevant data is retained in memory by the simulation, and
    provided methods exist (and are implemented) to extract the information required for probing.
    In practice, this can be achieved for imported components, where the spatial distribution of
    the source or medium is defined through a set of particles or cells. Because these quantities
    are probed directly from the original (imported) data structures, there are no (additional)
    grid effects.

    Forms describe how a given quantity should be probed. Generic forms (inheriting the GenericForm
    class) can be associated with any probe, including both spatial grid and input model probes.
    Examples include planar and linear cuts through the spatial domain, sampling at a
    user-specified list of positions, and parallel or all-sky projections. Spatial grid forms
    (inheriting the SpatialGridForm class) rely on the fact that the probed quantity is discretized
    on the simulation's spatial grid, and thus can be associated only with spatial grid probes. The
    prime example is the per-cell form.

    <b>The probe-form bridge</b>

    As can be deduced from the description above, the various types of probes and forms offer or
    require different types of primitives for retrieving probed quantity values. The
    ProbeFormBridge class provides a bridge between these differing capabilities and requirements.
    When a probe wants to output a file, it constructs a bridge instance and calls one of the
    bridge's writeQuantity() functions, passing it the appropriate information and local callback
    functions for retrieving probed quantity values on request. The bridge in turn calls the form's
    writeQuantity(bridge) function. Each Form subclass offers its own implementation of this
    function depending on the form's type of output. The form can get all required information
    through the bridge, including probed quantity values. The bridge passes the latter requests to
    the probe's call back functions, translating between request types as needed. For example, in
    the case of a spatial grid probe, a request by the form for a value at a given position is
    translated to a request for a value in the cell containing the position.

    One aspect that needs special care is the conversion from internal units (always SI) to output
    units (defined by the unit system configured by the user). For each probe-form combination, the
    bridge handles two different quantity types: the 'straight' quantity, probed at a given
    position, and the 'projected' quantity, probed along a given path. The projected quantity can
    be averaged or accumulated over the path. In the first case, the project and straight
    quantities have the same units; in the second case they don't.

    The bridge always supplies output units to the form, simplifying the implementation of a Form
    subclass. On the probe side, things are more complicated. Probes that handle scalar or vector
    quantities must specify the type of quantity to the bridge and supply any quantity values in
    internal units. The bridge will take care of the unit conversions. On the other hand, probes
    handling compound quantities must always convert quantity values to output units before
    supplying them to the bridge. This is so because the conversion factor may differ between
    components in a complicated way (for example, spectral quantities at different wavelengths).

    There are several variations of the writeQuantity() function tailored to scalar, vector and
    compound quantities, and catering to spatial grid or input model probes. Because these
    functions often expect the same arguments in different combinations, all arguments are
    described here.

    - \em fileid, \em projectedFileid: a string identifying the output file for the quantity being
    probed, straight and projected, among the possibly multiple files produced by this probe (e.g.
    "v", "1_rho", or "dust_T").

    - \em quantity, \em projectedQuantity: the quantity names (as defined in the SKIRT unit system)
    of the quantity being probed, straight and projected. These strings differ if the projected
    quantity is accumulated over distance and they are equal if the projected quantity is a
    weighted average. If these arguments are specified, the bridge performs all unit conversions.

    - \em unit, \em projectedUnit: the unit strings (as defined in the SKIRT unit system) of the
    quantity being probed, straight and projected. If these arguments are specified, the bridge
    does not perform any unit conversions and the probe is expected to provide values in output
    units.

    - \em description, \em projectedDescription: the user-oriented descriptions corresponding to
    the straight and projected quantities. These strings may differ even if the quantity names or
    unit strings are the same, for example, to indicate the weighting mechanism.

    - \em axis: for compound values, the values on the compound axis (e.g., wavelengths), in output
    units.

    - \em axisUnit: for compound values, the unit string for the compound axis values.

    - \em addColumnDefinitions: the callback function that will be called to add the appropriate
    column definitions if the form creates a text column output file.

    - \em valueInCell: the callback function that will be used to retrieve values of the quantity
    being probed in a spatial cell with given index. Used for spatial grid probes.

    - \em weightInCell: the callback function that will be used to retrieve the weight of the
    quantity being probed in a spatial cell with given index, in case averaging is required. Used
    for spatial grid probes.

    - \em valueAtPosition: the callback function that will be used to retrieve values of the
    quantity being probed at the specified position. Used for input model probes.

    - \em valueAlongPath: the callback function that will be used to retrieve values of the
    quantity being probed along the path starting at the specified position in the specified
    direction. Used for input model probes.

    - \em snapshots: a list of pointers to the snapshots representing the set of imported source
    or medium components being probed. Used for input model probes that probe imported source or
    medium components.

    - \em valueInEntity: the callback function that will be used to retrieve values of the quantity
    being probed in a snapshot entity with given index. Used for input model probes that probe an
    imported source or medium component.

    - \em weightInEntity: the callback function that will be used to retrieve the weight of the
    quantity being probed in a snapshot entity with given index, in case averaging is required.
    Used for input model probes that probe an imported source or medium component.

    */
class ProbeFormBridge
{
public:
    //======== Construction: for use by all probe types  =======

    /** The sole constructor accepts and retains pointers to the probe and form to be associated
        with this bridge. */
    ProbeFormBridge(const Probe* probe, const Form* form);

    //======== Callback type declarations  =======

    /** This is the type declaration for the callback function provided by the probe to add column
        definitions for a straight compound quantity to the specified text output file. */
    using AddColumnDefinitions = std::function<void(TextOutFile& outfile)>;

    /** This is the type declaration for the callback function provided by the spatial grid probe
        to retrieve the scalar value (in internal units) of the quantity being probed in the
        spatial cell with index \f$m\f$. */
    using ScalarValueInCell = std::function<double(int m)>;

    /** This is the type declaration for the callback function provided by the spatial grid probe
        to retrieve the vector value (in internal units) of the quantity being probed in the
        spatial cell with index \f$m\f$. */
    using VectorValueInCell = std::function<Vec(int m)>;

    /** This is the type declaration for the callback function provided by the spatial grid probe
        to retrieve the compound value (in output units) of the quantity being probed in the
        spatial cell with index \f$m\f$. The returned array must have the same number of elements
        as the \em axis array passed to the writeQuantity() function. */
    using CompoundValueInCell = std::function<Array(int m)>;

    /** This is the type declaration for the callback function provided by the spatial grid probe
        to retrieve the weight (in arbitrary units) of the quantity being probed in the spatial
        cell with index \f$m\f$. */
    using WeightInCell = std::function<double(int m)>;

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the scalar value (in internal units) of the quantity being probed at the specified
        position. */
    using ScalarValueAtPosition = std::function<double(Position bfr)>;

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the vector value (in internal units) of the quantity being probed at the specified
        position. */
    using VectorValueAtPosition = std::function<Vec(Position bfr)>;

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the compound value (in output units) of the quantity being probed at the specified
        position. The returned array must have the same number of elements as the \em axis array
        passed to the writeQuantity() function. */
    using CompoundValueAtPosition = std::function<Array(Position bfr)>;

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the scalar value (in internal units) of the quantity being probed along the path
        starting at the specified position in the specified direction. */
    using ScalarValueAlongPath = std::function<double(Position bfr, Direction bfk)>;

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the vector value (in internal units) of the quantity being probed along the path
        starting at the specified position in the specified direction. */
    using VectorValueAlongPath = std::function<Vec(Position bfr, Direction bfk)>;

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the compound value (in output units) of the quantity being probed along the path
        starting at the specified position in the specified direction. The returned array must have
        the same number of elements as the \em axis array passed to the writeQuantity() function.
        */
    using CompoundValueAlongPath = std::function<Array(Position bfr, Direction bfk)>;

    /** This is the type declaration for the callback function provided by the imported model probe
        to retrieve the scalar value (in internal units) of the quantity being probed in the
        snapshot entity with index \f$m\f$. */
    using ScalarValueInEntity = std::function<double(const Snapshot* snapshot, int m)>;

    /** This is the type declaration for the callback function provided by the imported model probe
        to retrieve the vector value (in internal units) of the quantity being probed in the
        snapshot entity with index \f$m\f$. */
    using VectorValueInEntity = std::function<Vec(const Snapshot* snapshot, int m)>;

    /** This is the type declaration for the callback function provided by the imported model probe
        to retrieve the weight (in arbitrary units) of the quantity being probed in the snapshot
        entity with index \f$m\f$. */
    using WeightInEntity = std::function<double(const Snapshot* snapshot, int m)>;

    //======== Writing: for use by spatial grid probes  =======

    /** This function causes the form associated with this bridge to output a file for a scalar
        quantity that needs to be accumulated along a path, according to the provided information.
        It should be called only from spatial grid probes. Refer to the class header for more
        information on the arguments. */
    void writeQuantity(string fileid, string projectedFileid, string quantity, string projectedQuantity,
                       string description, string projectedDescription, ScalarValueInCell valueInCell);

    /** This function causes the form associated with this bridge to output a file for a scalar
        quantity that needs to be averaged along a path, according to the provided information. It
        should be called only from spatial grid probes. Refer to the class header for more
        information on the arguments. */
    void writeQuantity(string fileid, string quantity, string description, string projectedDescription,
                       ScalarValueInCell valueInCell, WeightInCell weightInCell);

    /** This function causes the form associated with this bridge to output a file for a vector
        quantity (which is always averaged along a path) according to the provided information. It
        should be called only from spatial grid probes. Refer to the class header for more
        information on the arguments. */
    void writeQuantity(string fileid, string quantity, string description, string projectedDescription,
                       VectorValueInCell valueInCell, WeightInCell weightInCell);

    /** This function causes the form associated with this bridge to output a file for a compound
        quantity that needs to be accumulated along a path according to the provided information.
        It should be called only from spatial grid probes. Refer to the class header for more
        information on the arguments. */
    void writeQuantity(string fileid, string projectedFileid, string quantity, string projectedQuantity,
                       string description, string projectedDescription, const Array& axis, string axisUnit,
                       AddColumnDefinitions addColumnDefinitions, CompoundValueInCell valueInCell);

    /** This function causes the form associated with this bridge to output a file for a compound
        quantity that needs to be averaged along a path according to the provided information. It
        should be called only from spatial grid probes. Refer to the class header for more
        information on the arguments. */
    void writeQuantity(string fileid, string unit, string description, string projectedDescription, const Array& axis,
                       string axisUnit, AddColumnDefinitions addColumnDefinitions, CompoundValueInCell valueInCell,
                       WeightInCell weightInCell);

    //======== Writing: for use by input model probes  =======

    /** This function causes the form associated with this bridge to output a file for a scalar
        quantity according to provided the information. It should be called only from input model
        probes. Refer to the class header for more information on the arguments. */
    void writeQuantity(string fileid, string projectedFileid, string quantity, string projectedQuantity,
                       string description, string projectedDescription, ScalarValueAtPosition valueAtPosition,
                       ScalarValueAlongPath valueAlongPath);

    /** This function causes the form associated with this bridge to output a file for a vector
        quantity according to provided the information. It should be called only from input model
        probes. Refer to the class header for more information on the arguments. */
    void writeQuantity(string fileid, string projectedFileid, string quantity, string projectedQuantity,
                       string description, string projectedDescription, VectorValueAtPosition valueAtPosition,
                       VectorValueAlongPath valueAlongPath);

    /** This function causes the form associated with this bridge to output a file for a compound
        quantity according to provided the information. It should be called only from input model
        probes. Refer to the class header for more information on the arguments. */
    void writeQuantity(string fileid, string projectedFileid, string unit, string projectedUnit, string description,
                       string projectedDescription, const Array& axis, string axisUnit,
                       AddColumnDefinitions addColumnDefinitions, CompoundValueAtPosition valueAtPosition,
                       CompoundValueAlongPath valueAlongPath);

    //======== Writing: for use by imported source or media component probes  =======

    /** This function causes the form associated with this bridge to output a file for a scalar
        quantity that needs to be accumulated along a path, according to the provided information.
        It should be called only from input model probes for imported source or media components.
        Refer to the class header for more information on the arguments. */
    void writeQuantity(string fileid, string projectedFileid, string quantity, string projectedQuantity,
                       string description, string projectedDescription, const vector<const Snapshot*>& snapshots,
                       ScalarValueInEntity valueInEntity);

    /** This function causes the form associated with this bridge to output a file for a scalar
        quantity that needs to be averaged along a path, according to the provided information. It
        should be called only from input model probes for imported source or media components.
        Refer to the class header for more information on the arguments. */
    void writeQuantity(string fileid, string quantity, string description, string projectedDescription,
                       const vector<const Snapshot*>& snapshots, ScalarValueInEntity valueInEntity,
                       WeightInEntity weightInEntity);

    /** This function causes the form associated with this bridge to output a file for a vector
        quantity that needs to be averaged along a path, according to the provided information. It
        should be called only from input model probes for imported source or media components.
        Refer to the class header for more information on the arguments. */
    void writeQuantity(string fileid, string quantity, string description, string projectedDescription,
                       const vector<const Snapshot*>& snapshots, VectorValueInEntity valueInEntity,
                       WeightInEntity weightInEntity);

    //======== Querying: for use by all form types  =======

    /** This function returns the probe associated with this bridge upon construction. */
    const SimulationItem* probe() const;

    /** This function returns the simulation's spatial grid, or the null pointer if the simulation
        does not include any media. */
    const SpatialGrid* grid() const;

    /** This function returns the simulation's unit system. */
    const Units* units() const;

    /** This function returns the file name prefix formed by combining the probe name with the file
        identifier for the straight quantity set by the most recent call to writeQuantity(). */
    string prefix() const;

    /** This function returns the file name prefix formed by combining the probe name with the file
        identifier for the projected quantity set by the most recent call to writeQuantity(). */
    string projectedPrefix() const;

    /** This function returns the unit string of the straight quantity, as set by the most recent
        call to writeQuantity(). */
    string unit() const;

    /** This function returns the unit string of the projected quantity, as set by the most recent
        call to writeQuantity(). */
    string projectedUnit() const;

    /** This function returns the user description of the straight quantity, as set by the most
        recent call to writeQuantity(). */
    string description() const;

    /** This function returns the user description of the projected quantity, as set by the most
        recent call to writeQuantity(). */
    string projectedDescription() const;

    /** This function returns an array with the axis values for the quantity being probed, as set
        by the most recent call to writeQuantity(). For a scalar quantity, this returns the empty
        array (as opposed to an array with a single value). For a vector quantity, this returns an
        array containing {1.,2.,3.}. For a compound quantity, this returns the axis values set by
        the call to writeQuantity(). In the latter two cases, the number of elements in the
        returned array equals the return value of numValues(). */
    Array axis() const;

    /** This function returns the unit string of the axis values for the quantity, as set by the
        most recent call to writeQuantity(). For a scalar or vector quantity, this returns "1". For
        a compound quantity, this returns the axis unit set by the call to writeQuantity(). */
    string axisUnit() const;

    /** This function returns the number of values in the quantity being probed, as set by the most
        recent call to writeQuantity(). For a scalar quantity, this returns 1. For a vector
        quantity, this returns 3. For a compound quantity, this returns the size of the axis array
        set by the call to writeQuantity(). */
    int numValues() const;

    /** This function returns if the quantity being probed, as set by the most recent call to
        writeQuantity(), is a vector quantity, and false if it is a scalar or compound quantity. */
    bool isVector() const;

    /** This function adds the column definitions appropriate for the straight quantity being
        probed to the specified text output file. It should not be called for projected quantities.
        For scalar and vector quantities, the function adds the appropriate column definitions
        based on the description and unit information provided to the most recent call to
        writeQuantity(). For compound quantities, it invokes the \em addColumnDefinitions callback
        set by the most recent call to writeQuantity(). */
    void addColumnDefinitions(TextOutFile& outfile) const;

    //======== Querying: for use by spatial grid-specific forms  =======

    /** This function retrieves the values (in output units) of the quantity being probed in the
        spatial cell with index \f$m\f$. The provided array must have \em numValues elements that
        will be overwritten by this function.

        This function should be called only from within the writeQuantity() function of a spatial
        grid-specific form. It operates by invoking the \em valuesInCell callback set by the
        writeQuantity() function of this bridge that triggered the current file output. */
    void valuesInCell(int m, Array& values) const;

    //======== Querying: for use by generic forms  =======

    /** This function retrieves the values (in output units) of the quantity being probed at the
        specified position. The provided array must have \em numValues elements that will be
        overwritten by this function.

        This function should be called only from within the writeQuantity() function of a generic
        form. Depending on the type of the associated probe, it operates by invoking either the \em
        valuesAtPosition or \em valuesInCell callback set by the writeQuantity() function of this
        bridge that triggered the current file output. In the latter case, the function accesses
        the simulation's spatial grid to retrieve the index of the cell containing the specified
        position. */
    void valuesAtPosition(Position bfr, Array& values) const;

    /** This function retrieves the values (in output units) of the quantity being probed along the
        path starting at the specified position in the specified direction. The provided array must
        have \em numValues elements that will be overwritten by this function.

        This function should be called only from within the writeQuantity() function of a generic
        form. Its operation depends on the type of the associated probe, the quantity properties
        set by the most recent call to writeQuantity(), and the callbacks set by the
        writeQuantity() function of this bridge that triggered the current file output.

        For an input model probe, the function simply invokes the \em valuesAlongPath callback. For
        a spatial grid probe, the function accesses the simulation's spatial grid to determine the
        list of path segments (traversing each cell) corresponding to the specified path, and it
        calls the \em valuesInCell callback for each cell in the list. If the straight and
        projected quantity names differ, the function calculates the accumulated quantity
        \f$\sum_m Q_m \Delta s_m\f$. If the quantity names are the same, it calculates the
        weighted average \f$\sum_m Q_m w_m \Delta s_m / \sum_m w_m \Delta s_m\f$, where the
        weights are retrieved by calling the \em weightInCell callback for each cell. */
    void valuesAlongPath(Position bfr, Direction bfk, Array& values) const;

    //======== Data members: private  =======

private:
    // class instances -- set by constructor
    const Probe* _probe{nullptr};
    const Form* _form{nullptr};
    const SpatialGrid* _grid{nullptr};
    const Units* _units{nullptr};

    // probe/quantity type
    enum class Type {
        GridScalarAccumulated,
        GridScalarAveraged,
        GridVectorAveraged,
        GridCompoundAccumulated,
        GridCompoundAveraged,
        InputScalar,
        InputVector,
        InputCompound
    };
    Type _type{Type::GridScalarAccumulated};

    // quantity info -- set by writeQuantity()
    string _fileid;
    string _projectedFileid;
    string _unit;
    string _projectedUnit;
    double _unitFactor{1.};
    double _projectedUnitFactor{1.};
    string _description;
    string _projectedDescription;
    Array _axis;
    string _axisUnit;
    int _numValues{0};

    // call-back functions -- set by writeQuantity()
    AddColumnDefinitions _addColumnDefinitions;
    ScalarValueInCell _scalarValueInCell;
    VectorValueInCell _vectorValueInCell;
    CompoundValueInCell _compoundValueInCell;
    WeightInCell _weightInCell;
    ScalarValueAtPosition _scalarValueAtPosition;
    VectorValueAtPosition _vectorValueAtPosition;
    CompoundValueAtPosition _compoundValueAtPosition;
    ScalarValueAlongPath _scalarValueAlongPath;
    VectorValueAlongPath _vectorValueAlongPath;
    CompoundValueAlongPath _compoundValueAlongPath;
};

//////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBEFORMBRIDGE_HPP
#define PROBEFORMBRIDGE_HPP

#include "Array.hpp"
#include <functional>
class Direction;
class Form;
class Position;
class Probe;
class TextOutFile;

//////////////////////////////////////////////////////////////////////

/** ProbeFormBridge is a helper class that forms a bridge between Probe and Form instances for
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
    In practice, this can be achieved only in specific cases, e.g. for the set of smoothed
    particles defining the spatial distribution of a primary source. Because these quantities are
    probed directly from the original (imported) data structures, there are no (additional) grid
    effects.

    Forms describe how a given quantity should be probed. Generic forms (inheriting the GenericForm
    class) can be associated with any probe, including both spatial grid and input model probes.
    Examples include planar and linear cuts through the spatial domain, sampling at a
    user-specified list of positions, and parallel or all-sky projections. Spatial grid forms
    (inheriting the SpatialGridForm class) rely on the fact that the probed quantity is discretized
    on the simulation's spatial grid, and thus can be associated only with spatial grid probes. The
    prime (and perhaps only) example is the per-cell form.

    <b>The probe-form bridge</b>

    As can be deduced from the description above, the various types of probes and forms offer or
    require different types of primitives for retrieving probed quantity values. The
    ProbeFormBridge class provides a bridge between these differing capabilities and requirements.
    When a probe wants to output a file, it constructs and configures a bridge instance and calls
    one of the bridge's writeFile() functions, passing it the appropriate local callback functions
    for retrieving probed quantity values on request. The bridge in turn calls the form's
    writeFile(bridge) function. Each Form subclass offers its own implementation of this function
    depending on the form's type of output. The form can get all required information through the
    bridge, including probed quantity values. The bridge passes the latter requests to the probe's
    call back functions, translating between request types as needed. For example, in the case of a
    spatial grid probe, a request by the form for a value at a given position is translated to a
    request for a value in the cell containing the position. */
class ProbeFormBridge
{
public:
    //======== Construction and setup: for use by all probe types  =======

    /** The sole constructor accepts and retains pointers to the probe and form to be associated
        with this bridge. */
    ProbeFormBridge(const Probe* probe, const Form* form);

    /** This is the type declaration for the callback function provided by the probe to add the
        appropriate column definitions the specified text output file. */
    using AddColumnHeaders = std::function<void(TextOutFile* outfile)>;

    /** This function sets the properties of the quantity being probed. The probe should call
        setQuantity() before invoking writeFile(). In most cases, a single call is sufficient, even
        if muliple files are being written. However, multiple calls are allowed. For example, a
        given probe may produce mass densities for dust and number densities for gas media.

        The specified properties include:

        - The number of component values in the quantity (e.g., 1 for density, 3 for velocity,
        \f$N_\lambda\f$ for the radiation field).

        - The quantity names (as defined for the SKIRT unit system) of the quantity being probed,
        straight (as returned by the valuesAtPosition and \em valuesInCell callbacks) and projected
        (as returned by the \em valuesAlongPath callback). These strings differ if the projected
        quantity is accumulated over distance and they are equal if the projected quantity is a
        weighted average.

        - The user-oriented descriptions corresponding to the straight and projected quantities.
        These strings may differ even if the quantity names are the same, for example, to indicate
        the weighting mechanism.

        - The callback function that will be called to add the appropriate column definitions if
        the form creates a text output file in response to a writeFile() invocation. The callback
        must add exactly \em numValues column definitions. */
    void setQuantity(int numValues, string quantity, string projectedQuantity, string description,
                     string projectedDescription, AddColumnHeaders addColumnHeaders);

    //======== Writing: for use by spatial grid probes  =======

    /** This is the type declaration for the callback function provided by the spatial grid probe
        to retrieve the values of the quantity being probed in the spatial cell with index \f$m\f$.
        The provided array will have \em numValues elements with undefined values that need to be
        overwritten by the callback. */
    using ValuesInCell = std::function<void(int m, Array& values)>;

    /** This is the type declaration for the callback function provided by the spatial grid probe
        to retrieve the weight of the quantity being probed in the spatial cell with index \f$m\f$.
        The callback is used only in case the straight and projected quantity names differ. */
    using WeightInCell = std::function<double(int m)>;

    /** This function causes the form associated with this bridge to output a file for the quantity
        being probed according to the most recent call to setQuantity(). It should be called only
        from spatial grid probes.

        The specified arguments include:

        - a string identifying the file among the possibly multiple files produced by this probe
        (e.g. "0", "1", or "dust"), or the empty string if such identification is not needed.

        - the callback function that will be used to retrieve values of the quantity being probed
        in the spatial cell with index \f$m\f$. The provided array will have \em numValues elements
        with undefined values that need to be overwritten by the callback.

        - the callback function that will be used to retrieve the weight of the quantity being
        probed in the spatial cell with index \f$m\f$. This callback is used only in case the
        straight and projected quantity names are equal. If the quantity names differ, the argument
        may be omitted. */
    void writeFile(string fileid, ValuesInCell valuesInCell, WeightInCell weightInCell = nullptr);

    //======== Writing: for use by input model probes  =======

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the values of the quantity being probed at the specified position. The provided
        array will have \em numValues elements with undefined values that need to be overwritten by
        the callback. */
    using ValuesAtPosition = std::function<void(Position bfr, Array& values)>;

    /** This is the type declaration for the callback function provided by the input model probe to
        retrieve the values of the quantity being probed along the path starting at the specified
        position in the specified direction. The provided array will have \em numValues elements
        with undefined values that need to be overwritten by the callback. */
    using ValuesAlongPath = std::function<double(Position bfr, Direction bfk, Array& values)>;

    /** This function causes the form associated with this bridge to output a file for the quantity
        being probed according to the most recent call to setQuantity(). It should be called only
        from input model probes.

        The specified arguments include:

        - a string identifying the file among the possibly multiple files produced by this probe
        (e.g. "0", "1", or "dust"), or the empty string if such identification is not needed.

        - the callback function that will be used to retrieve values of the quantity being probed
        at the specified position. The provided array will have \em numValues elements with
        undefined values that need to be overwritten by the callback.

        - the callback function that will be used to retrieve values of the quantity being probed
        along the path starting at the specified position in the specified direction. The provided
        array will have \em numValues elements with undefined values that need to be overwritten by
        the callback. */
    void writeFile(string fileid, ValuesAtPosition valuesAtPosition, ValuesAlongPath valuesAlongPath);

    //======== Querying: for use by all form types  =======

    /** This function returns the probe associated with this bridge upon construction. */
    Probe* probe() const;

    /** This function returns the number of values in the quantity being probed, as set by the most
        recent call to setQuantity(). */
    int numValues() const;

    /** This function returns the name of the straight quantity, as set by the most recent call to
        setQuantity(). */
    string quantity() const;

    /** This function returns the name of the projected quantity, as set by the most recent call to
        setQuantity(). */
    string projectedQuantity() const;

    /** This function returns the user description of the straight quantity, as set by the most
        recent call to setQuantity(). */
    string description() const;

    /** This function returns the user description of the projected quantity, as set by the most
        recent call to setQuantity(). */
    string projectedDescription() const;

    /** This function adds the column definitions appropriate for the quantity being probed to the
        specified text output file. It does so by invoking the \em addColumnHeaders callback set by
        the most recent call to setQuantity(). */
    void addColumnHeaders(TextOutFile* outfile) const;

    //======== Querying: for use by spatial grid-specific forms  =======

    /** This function retrieves the values of the quantity being probed in the spatial cell with
        index \f$m\f$. The provided array must have \em numValues elements with undefined values
        that will be overwritten by this function.

        This function should be called only from within the writeFile() function of a spatial
        grid-specific form. It operates by invoking the \em valuesInCell callback set by the
        writeFile() function of this bridge that triggered the current file output. */
    void valuesInCell(int m, Array& values) const;

    //======== Querying: for use by generic forms  =======

    /** This function retrieves the values of the quantity being probed at the specified position.
        The provided array must have \em numValues elements with undefined values that will be
        overwritten by this function.

        This function should be called only from within the writeFile() function of a generic form.
        Depending on the type of the associated probe, it operates by invoking either the \em
        valuesAtPosition or \em valuesInCell callback set by the writeFile() function of this
        bridge that triggered the current file output. In the latter case, the function accesses
        the simulation's spatial grid to retrieve the index of the cell containing the specified
        position. */
    void valuesAtPosition(Position bfr, Array& values) const;

    /** This function retrieves the values of the quantity being probed along the path starting at
        the specified position in the specified direction. The provided array must have \em
        numValues elements with undefined values that will be overwritten by this function.

        This function should be called only from within the writeFile() function of a generic form.
        Its operation depends on the type of the associated probe, the quantity properties set by
        the most recent call to setQuantity(), and the callbacks set by the writeFile() function of
        this bridge that triggered the current file output.

        For an input model probe, the function simply invokes the \em valuesAlongPath callback. For
        a spatial grid probe, the function accesses the simulation's spatial grid to determine the
        list of path segments (traversing each cell) corresponding to the specified path, and it
        calls the \em valuesInCell callback for each cell in the list. If the straight and
        projected quantity names differ, the function calculates the accumulated quantity
        \f$\Sigma_m Q_m \Delta s_m\f$. If the quantity names are the same, it calculates the
        weigthed average \f$\Sigma_m Q_m w_m \Delta s_m / \Sigma_m w_m \Delta s_m\f$, where the
        weigths are retrieved by calling the \em weightInCell callback for each cell. */
    void valuesAlongPath(Position bfr, Direction bfk, Array& values) const;

    //======== Data members: private  =======

private:
    // set by constructor
    const Probe* _probe{nullptr};
    const Form* _form{nullptr};

    // set by setQuantity()
    int _numValues{0};
    string _quantity;
    string _projectedQuantity;
    string _description;
    string _projectedDescription;
    AddColumnHeaders _addColumnHeaders;
};

//////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOREDTABLE_HPP
#define STOREDTABLE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include "NR.hpp"
#include "Range.hpp"
#include "StoredTableImpl.hpp"
#include <array>

////////////////////////////////////////////////////////////////////

/** An instance of the StoredTable<N> class template provides access to the contents of a
    particular resource file in the SKIRT stored table format (i.e. a "stored table").

    Stored table file format
    ------------------------

    A stored table includes the names of the quantities on the axes (e.g. wavelength and grain
    size) and those being tabulated (e.g. absorption and scattering efficiencies), in addition to
    the grid points for each axis and the tabulated data values. All grid points and values are
    stored as binary data in the form of 64-bit floating-point numbers, and are always given in SI
    units. The format is designed so that it is easy to calculate the offset, relative to the start
    of the file, to any particular data value. More specifically, a stored table file is
    essentially a sequence of 8-byte data items. A data item can have one of three types:
        - string: 1 to 8 printable and non-whitespace 7-bit ASCII characters, padded with spaces
          to fill 8 bytes if needed;
        - unsigned integer: 64-bit integer in little-endian byte order;
        - floating point: 64-bit double (IEEE 754) in little-endian byte order.

    The overall layout is as follows:
        - SKIRT name/version tag
        - Endianness tag
        - numAxes
        - axisName (x numAxes)
        - axisUnit (x numAxes)
        - axisScale (x numAxes)
        - [ numPoints  axisPoint (x numPoints) ] (x numAxes)
        - numQuantities
        - quantityName (x numQuantities)
        - quantityUnit (x numQuantities)
        - quantityScale (x numQuantities)
        - value (x numQuantities x numPoints1 x ... x numPointsN)
        - end-of-file tag

    The values are ordered so that the quantity values for a particular point are next to each
    other, the first axis index varies most rapidly, and the last axis index varies least rapidly.

    The StoredTable<N> class template
    ---------------------------------

    The StoredTable<N> template parameter \em N specifies the number of axes in the stored table,
    and thus the number of axis values needed to retrieve a tabulated quantity from the table. Each
    StoredTable<N> instance represents a single tabulated quantity. Accessing the full contents of
    a stored table resource file with multiple tabulated quantities requires a seperate
    StoredTable<N> instance for each of those quantities.

    The default constructor creates an invalid stored table instance. The alternate constructor and
    the open() function associate a particular stored table resource file with the stored table
    instance. The number of axes in this stored table resource file must match the template
    parameter \em N. Also, the axis names and the corresponding units in the file, and one of the
    tabulated quantity names and its corresponding unit in the file, must match the information
    passed to the alternate constructor or the open() function. The destructor automatically
    releases the file association and any related resources.

    The parenthesis operator returns the quantity represented by the StoredTable<N> for the \em N
    specified axis values, interpolated from the tabulated values. Other functions offer access
    for specific purposes, such as constructing a cumulative distribution function along one axis,
    given values for the other axes.

    Implementation and performance
    ------------------------------

    A StoredTable<N> instance acquires a read-only memory map on the associated stored table
    resource file as opposed to actually reading the file contents into memory through regular file
    I/O operations. This has some important, mostly positive, consequences.

    Acquiring a memory map establishes a mapping between "pages" of system-defined size in the
    logical address space of a process and the contents of the "backing file", in this case the
    stored table resource file. This operation is simple and thus very fast. From then on, the
    operating system automatically loads pages from the backing file into physical memory as they
    become needed because the program addresses an item in the logical memory range of the page.
    Conversely, the operating system automatically removes pages from physical memory if available
    memory becomes tight. In effect, the operating system automatically manages a high-performance
    caching mechanism on stored tables.

    Three important use cases benefit greatly from this mechanism. Firstly, a large resource file
    can be left associated with a StoredTable<N> instance for the duration of the program, even if
    it is used only sporadically. When memory is tight, infrequently used portions of the data will
    automatically be removed from memory and reloaded later if needed. Secondly, there is little
    overhead in constructing a StoredTable<N> instance (and possibly destroying it shortly
    thereafter) even when the program needs only a small portion of the file contents. And thirdly,
    because all StoredTable<N> instances associated with a given stored table resource file share
    the same memory map on that file, using a seperate instance for each quantity in the stored
    table incurs very little overhead.

    Moreover, most operating systems share memory maps between processes. For a parallel program
    using MPI, this means that all processes running on the same compute node share a single memory
    copy of the resources they employ. Also, most operating systems keep the memory map caches
    alive between consecutive invocations of a program (assuming memory is available), increasing
    performance when, for example, interactively testing the program.

    On the downside, a program requesting a huge chunk of data from a large stored table in a
    serial fashion would run faster using regular file I/O, because the separate page loads take
    more time than sequentially reading data in bulk. More importantly, performance usually
    degrades rapidly (to the point where the program no longer performs any useful work) when the
    system is forced to constantly remove and reload pages because there is not enough memory to
    hold the data needed for a particular phase in the program. And finally, the run-time
    performance of a program becomes somewhat unpredicable because the speed of accessing resources
    depends heavily on the previous state of the operating system caches.
*/
template<size_t N> class StoredTable
{
    static_assert(N >= 1, "StoredTable number of axes must be at least 1");

    // ================== Constructing ==================

public:
    /** The default constructor constructs an invalid stored table instance. The user must call the
        open() function to associate the stored table instance with a particular stored table
        resource file. Calling any of the other functions before calling open() results in
        undefined behavior (usually a crash). */
    StoredTable() {}

    /** This alternate constructor constructs a stored table instance and immediately associates a
        given stored table resource file with it by calling the open() function. Refer to the
        open() function for a description of the arguments and of its operation. */
    StoredTable(const SimulationItem* item, string filename, string axes, string quantity, bool clampFirstAxis = true,
                bool resource = true)
    {
        open(item, filename, axes, quantity, clampFirstAxis, resource);
    }

    /** This function associates a given stored table resource or input file with the stored table
        instance. If such an association already exists, this function throws a fatal error.
        Conversely, calling any of the other functions before an association exists results in
        undefined behavior (usually a crash).

        The \em item argument specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve an appropriate logger.

        If the \em resource flag is true (the default value), the \em filename argument specifies
        the filename of the built-in resource, without any directory segments. If the \em resource
        flag is false, the \em filename argument specifies the file path of the stored table input
        file relative to the input path of the simulation. In both cases, the file must have the
        ".stab" filename extension, which will be added to the specified filename if needed.

        The \em filename argument specifies the filename of the resource, without any directory
        segments. The resource file must have the ".stab" filename extension, which will be added
        to the specified filename if needed.

        First of all, the number of axes in this stored table resource file must match the template
        parameter \em N. Furthermore, the axes names in the resource file and the corresponding
        units for each axis must match the information specified in the \em axes argument. Finally,
        one of the tabulated quantity names in the resource file and its corresponding unit must
        match the information specified in the \em quantity argument. For a stored table resource
        file with multiple tabulated quantities, the \em quantity argument at the same time
        determines which of these quantities will be associated with the stored table instance.

        The string passed to the \em axes argument must have the syntax
        "name1(unit1),...,nameN(unitN)". In other words, a unit string between parenthesis follows
        each axis name, and the specifications for different axes are separated by a comma. For
        example, "lambda(m),a(m)". Whitespace is not allowed. The string passed to the \em quantity
        argument must have a similar syntax, for a single name/unit combination. Examples include
        "Llambda(W/m)", "Qabs(1)", and "h(J/m3)".

        By default, out-of-range axes values are clamped to the corresponding outer grid point. In
        other words, a "nearby" quantity value at some outer grid point is used for out-of-range
        axes values. This behavior can be turned off for the first axis (but not for the other
        axes) by passing a value of false to the optional \em clampFirstAxis flag. In that case,
        the quantity value for out-of-range first-axis values is considered to be zero.

        In summary, this function (1) locates the specified stored table resource file, (2)
        acquires a memory map on the file, (3) verifies that the stored table matches all
        requirements, and (4) stores relevant information in data members. If any of these steps
        fail, the function throws a fatal error. */
    void open(const SimulationItem* item, string filename, string axes, string quantity, bool clampFirstAxis = true,
              bool resource = true)
    {
        StoredTable_Impl::open(N, item, filename, resource, axes, quantity, _filePath, &_axBeg[0], &_qtyBeg, &_axLen[0],
                               &_qtyStep, &_axLog[0], &_qtyLog);
        _clamp = clampFirstAxis;
    }

    /** The destructor breaks the association with a stored table resource file established by the
        alternate constructor or the open() function, if there is any. In practice, this simply
        means releasing the memory map on the associated file. */
    ~StoredTable() { StoredTable_Impl::close(_filePath); }

    // ================== Accessing values ==================

public:
    /** This function returns the value of the quantity represented by this stored table for the
        specified axes values, interpolated over the grid points of the actual tabulated values in
        all dimensions. The function uses linear or logarithmic interpolation for the axes and
        quantity values according to the flags specified in the stored table. Out-of-range axes
        values are handled according to the policy set by the \em clampFirstAxis flag in the open()
        function. */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isFloatArgList<N, Values...>()>>
    double operator()(Values... values) const
    {
        // storage for each axis
        std::array<double, N> value = {{static_cast<double>(values)...}};
        std::array<size_t, N> i2;  // upper grid bin boundary index
        std::array<double, N> f;   // fraction of axis value in bin

        // precompute for each axis
        for (size_t k = 0; k != N; ++k)
        {
            // get the index of the upper border of the axis grid bin containing the specified axis value
            double x = value[k];
            size_t right = std::lower_bound(_axBeg[k], _axBeg[k] + _axLen[k], x) - _axBeg[k];

            // if the value is beyond the grid borders:
            //    - if we're not clamping, simply return zero
            //    - if we're clamping, adjust both the bin border and the value
            if (right == 0)
            {
                if (!_clamp && k == 0 && x != _axBeg[k][0]) return 0.;
                right++;
                x = _axBeg[k][0];
            }
            else if (right == _axLen[k])
            {
                if (!_clamp && k == 0) return 0.;
                right--;
                x = _axBeg[k][right];
            }
            i2[k] = right;

            // get the axis values at the grid borders
            double x1 = _axBeg[k][right - 1];
            double x2 = _axBeg[k][right];

            // if requested, compute logarithm of coordinate values
            if (_axLog[k])
            {
                x = log(x);
                x1 = log(x1);
                x2 = log(x2);
            }

            // calculate the fraction of the requested axis value in the bin
            f[k] = (x - x1) / (x2 - x1);
        }

        // storage for each term in the interpolation
        constexpr size_t numTerms = 1 << N;  // there are 2^N terms
        std::array<double, numTerms> ff;     // front factor
        std::array<double, numTerms> yy;     // tabulated value

        // determine front factor and tabulated value for each term of the interpolation
        std::array<size_t, N> indices;  // storage for indices of the current term
        for (size_t t = 0; t != numTerms; ++t)
        {
            // use the binary representation of the term index to determine left/right for each axis
            size_t term = t;  // temporary version of term index that will be bit-shifted
            double front = 1.;
            for (size_t k = 0; k != N; ++k)
            {
                size_t left = term & 1;  // lowest significant digit = 1 means lower border
                indices[k] = i2[k] - left;
                front *= left ? (1 - f[k]) : f[k];
                term >>= 1;
            }
            if (front)
            {
                ff[t] = front;
                yy[t] = valueAtIndices(indices);

                // if logarithmic interpolation of y value is requested and not all bordering values are positive,
                // we can't comply with the request, so return zero
                if (_qtyLog && yy[t] <= 0) return 0.;
            }
            else
            {
                ff[t] = 0.;
                yy[t] = 0.;
            }
        }

        // calculate sum
        double y = 0.;
        if (_qtyLog)
        {
            for (size_t t = 0; t != numTerms; ++t)
                if (ff[t]) y += ff[t] * log(yy[t]);
            y = exp(y);
        }
        else
        {
            for (size_t t = 0; t != numTerms; ++t) y += ff[t] * yy[t];
        }
        return y;
    }

    /** For a one-dimensional table only, this function returns the value of the quantity
        represented by the stored table for the specified axes value, interpolated over the grid
        points of the actual tabulated values. The function uses linear or logarithmic
        interpolation for the axis and quantity values according to the flags specified in the
        stored table. Out-of-range axes values are handled according to the policy set by the \em
        clampFirstAxis flag in the open() function. */
    template<typename Value, typename = std::enable_if_t<N == 1 && CompileTimeUtils::isFloatArgList<1, Value>()>>
    double operator[](Value value) const
    {
        return operator()(value);
    }

    /** This function returns a sequence of values of the quantity represented by this stored table
        (in the \em yv argument) corresponding to the specified sequence of first-axis values (the
        \em xv argument), using given fixed values for the other axes, if any (the arguments at the
        end of the list). The function behaves as if it would call the operator()() function for
        each element in the specified input sequence, but it is more efficient. */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isFloatArgList<N - 1, Values...>()>>
    void valueArray(Array& yv, const Array& xv, Values... values) const
    {
        // storage for each axis
        std::array<double, N> value = {{0., static_cast<double>(values)...}};
        std::array<size_t, N> i2;       // upper grid bin boundary index
        std::array<double, N> f;        // fraction of axis value in bin
        std::array<size_t, N> indices;  // indices of the current term

        // storage for each term in the interpolation
        constexpr size_t numTerms = 1 << N;  // there are 2^N terms
        std::array<double, numTerms> ff;     // front factor
        std::array<double, numTerms> yy;     // tabulated value

        // resize output array
        size_t n = xv.size();  // number of values to calculate
        yv.resize(n);

        // precompute grid index and fraction for all but the first axis
        for (size_t k = 1; k != N; ++k)
        {
            // get the index of the upper border of the axis grid bin containing the specified axis value
            double x = value[k];
            size_t right = std::lower_bound(_axBeg[k], _axBeg[k] + _axLen[k], x) - _axBeg[k];

            // if the value is beyond the grid borders, adjust both the bin border and the value
            if (right == 0)
            {
                right++;
                x = _axBeg[k][0];
            }
            else if (right == _axLen[k])
            {
                right--;
                x = _axBeg[k][right];
            }
            i2[k] = right;

            // get the axis values at the grid borders
            double x1 = _axBeg[k][right - 1];
            double x2 = _axBeg[k][right];

            // if requested, compute logarithm of coordinate values
            if (_axLog[k])
            {
                x = log(x);
                x1 = log(x1);
                x2 = log(x2);
            }

            // calculate the fraction of the requested axis value in the bin
            f[k] = (x - x1) / (x2 - x1);
        }

        // calculate each value in the array
        for (size_t i = 0; i != n; ++i)
        {
            value[0] = xv[i];

            // compute grid index and fraction for the first axis
            {
                // get the index of the upper border of the axis grid bin containing the specified axis value
                double x = value[0];
                size_t right = std::lower_bound(_axBeg[0], _axBeg[0] + _axLen[0], x) - _axBeg[0];

                // if the value is beyond the grid borders:
                //    - if we're not clamping, simply return zero
                //    - if we're clamping, adjust both the bin border and the value
                if (right == 0)
                {
                    if (!_clamp && x != _axBeg[0][0])
                    {
                        yv[i] = 0.;
                        continue;
                    };
                    right++;
                    x = _axBeg[0][0];
                }
                else if (right == _axLen[0])
                {
                    if (!_clamp)
                    {
                        yv[i] = 0.;
                        continue;
                    };
                    right--;
                    x = _axBeg[0][right];
                }
                i2[0] = right;

                // get the axis values at the grid borders
                double x1 = _axBeg[0][right - 1];
                double x2 = _axBeg[0][right];

                // if requested, compute logarithm of coordinate values
                if (_axLog[0])
                {
                    x = log(x);
                    x1 = log(x1);
                    x2 = log(x2);
                }

                // calculate the fraction of the requested axis value in the bin
                f[0] = (x - x1) / (x2 - x1);
            }

            // determine front factor and tabulated value for each term of the interpolation
            for (size_t t = 0; t != numTerms; ++t)
            {
                // use the binary representation of the term index to determine left/right for each axis
                size_t term = t;  // temporary version of term index that will be bit-shifted
                double front = 1.;
                for (size_t k = 0; k != N; ++k)
                {
                    size_t left = term & 1;  // lowest significant digit = 1 means lower border
                    indices[k] = i2[k] - left;
                    front *= left ? (1 - f[k]) : f[k];
                    term >>= 1;
                }
                if (front)
                {
                    ff[t] = front;
                    yy[t] = valueAtIndices(indices);

                    // if logarithmic interpolation of y value is requested and not all bordering values are positive,
                    // we can't comply with the request, so return zero
                    if (_qtyLog && yy[t] <= 0)
                    {
                        yv[i] = 0.;
                        continue;
                    };
                }
                else
                {
                    ff[t] = 0.;
                    yy[t] = 0.;
                }
            }

            // calculate sum
            double y = 0.;
            if (_qtyLog)
            {
                for (size_t t = 0; t != numTerms; ++t)
                    if (ff[t]) y += ff[t] * log(yy[t]);
                y = exp(y);
            }
            else
            {
                for (size_t t = 0; t != numTerms; ++t) y += ff[t] * yy[t];
            }
            yv[i] = y;
        }
    }

    // ------------------------------------------

    /** This function constructs both the normalized probability density function (pdf) and the
        corresponding normalized cumulative distribution function (cdf) for the tabulated quantity
        across a given range in the first axis (the \em xrange argument), using given fixed values
        for the other axes, if any (the arguments at the end of the list).

        The resulting first-axis grid is constructed into \em xv, the corresponding pdf into \em
        pv, and the corresponding cdf into \em Pv. In all cases, xv[0]=xmin, xv[n]=xmax, Pv[0]=0,
        and Pv[n]=1. The function returns the normalization factor, i.e. the value of Pv[n] before
        normalization.

        If any of the axes values, including the \em xrange values specifying the range for the
        first axis, are out of range of the internal grid, extra quantity values are fabricated
        according to the policy set by the \em clampFirstAxis flag in the open() function.

        If the flags specified in the stored table indicate that both the first axis and the
        quantity represented by the table should be interpolated logarithmically, it is assumed
        that the pdf behaves as a power-law between any two grid points, and the integration to
        determine the cdf is performed accordingly. In all other cases, piece-wise linear behavior
        is assumed and regular trapezium-rule integration is used. */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isFloatArgList<N - 1, Values...>()>>
    double cdf(Array& xv, Array& pv, Array& Pv, Range xrange, Values... values) const
    {
        // copy the relevant portion of the internal axis grid
        size_t minRight = std::upper_bound(_axBeg[0], _axBeg[0] + _axLen[0], xrange.min()) - _axBeg[0];
        size_t maxRight = std::lower_bound(_axBeg[0], _axBeg[0] + _axLen[0], xrange.max()) - _axBeg[0];
        xv.resize(2 + maxRight - minRight);  // number of internal grid points plus two external points
        size_t i = 0;                        // i = index in target array
        xv[i++] = xrange.min();              // j = index in internal grid array
        for (size_t j = minRight; j < maxRight;) xv[i++] = _axBeg[0][j++];
        xv[i++] = xrange.max();

        // interpolate or copy the corresponding probability density values
        valueArray(pv, xv, values...);

        // perform the rest of the operation in a non-templated function
        return NR::cdf2(_axLog[0] && _qtyLog, xv, pv, Pv);
    }

    // ------------------------------------------

    /** This function returns the range of the table axis indicated by the zero-based index in the
        template argument. */
    template<size_t axisIndex, typename = std::enable_if_t<axisIndex <= N>> Range axisRange() const
    {
        return Range(_axBeg[axisIndex][0], _axBeg[axisIndex][_axLen[axisIndex] - 1]);
    }

    // ================== Accessing the raw data in a 1D table ==================

public:
    /** This function is available only for one-dimensional tables. It returns the number of
        entries in the table, i.e. the number of grid points in the single axis and the number of
        corresponding quantity values. */
    size_t size() const
    {
        static_assert(N == 1, "This function is available only for one-dimensional tables");
        return _axLen[0];
    }

    /** This function is available only for one-dimensional tables. It returns a pointer to the
        first element in the axis data array. The number of elements in this array can be obtained
        through the size() function. */
    const double* axisData() const
    {
        static_assert(N == 1, "This function is available only for one-dimensional tables");
        return _axBeg[0];
    }

    /** This function is available only for one-dimensional tables. It returns a pointer to the
        first element in the quantity data array. The number of elements in this array can be
        obtained through the size() function. */
    const double* quantityData() const
    {
        static_assert(N == 1, "This function is available only for one-dimensional tables");
        return _qtyBeg;
    }

    // ================== Accessing the raw data ==================

private:
    /** This template function returns a copy of the value at the specified N indices. There is no
        range checking. Out-of-range index values cause unpredictable behavior. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isFloatArgList<N, Indices...>()>>
    double valueAtIndices(Indices... indices) const
    {
        return _qtyBeg[flattenedIndex(std::array<size_t, N>({{static_cast<size_t>(indices)...}}))];
    }

    /** This function returns a copy of the value at the specified N indices. There is no range
        checking. Out-of-range index values cause unpredictable behavior. */
    double valueAtIndices(const std::array<size_t, N>& indices) const { return _qtyBeg[flattenedIndex(indices)]; }

    /** This function returns the flattened index in the underlying data array for the specified N
        indices. */
    size_t flattenedIndex(const std::array<size_t, N>& indices) const
    {
        size_t result = indices[N - 1];
        for (size_t k = N - 2; k < N; --k) result = result * _axLen[k] + indices[k];
        return result * _qtyStep;
    }

    // ================== Data members ==================

private:
    string _filePath;                     // the canonical path to the associated stored table file
    std::array<const double*, N> _axBeg;  // pointer to first grid point for each axis
    const double* _qtyBeg;                // pointer to first quantity value
    std::array<size_t, N> _axLen;         // number of grid points for each axis
    size_t _qtyStep;                      // step size from one quantity value to the next (1=adjacent)
    std::array<bool, N> _axLog;           // interpolation type (true=log, false=linear) for each axis
    bool _qtyLog;                         // interpolation type (true=log, false=linear) for quantity
    bool _clamp;                          // value for out-of-range first-axis indices: true=clamped, false=zero
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ARRAYTABLE_HPP
#define ARRAYTABLE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include <array>
#include <functional>

////////////////////////////////////////////////////////////////////

/** The ArrayTable<N> class template implements multi-dimensional tables with special support for
    \em rows in the last table dimension. Specifically, a ArrayTable<N> instance holds an
    N-dimensional table of double values. The combination of the first N-1 indices addresses a row
    in an (N-1)-dimensional table of rows; the last index addresses a column within a row. The
    values are stored as a list of Array objects, allowing read and write access to individual
    values as well as to rows of values as a whole.

    It is possible and allowed to resize the rows individually through their reference. In that
    case the caller should ensure that all rows end up with the same size (the ArrayTable
    implementation does not check this requirement and does not rely on it). Note that the
    size(N-1) function returns the original row size as it has been specified in the constructor or
    in the most recent resize() invocation. The rowsize() function returns the size of the first
    row in the table, or zero if the table is empty. */
template<size_t N> class ArrayTable
{
    static_assert(N >= 2, "ArrayTable dimension must be at least 2");

    // ================== Constructing ==================

public:
    /** The default constructor constructs an empty table. */
    ArrayTable() { _sizes.fill(0); }

    /** This constructor constructs a table holding the specified number of items in each of the N
        dimensions. If the last dimension in the list is nonzero, the arrays holding the data are
        resized correspondingly and all values are set to zero. If the last dimension is zero, the
        arrays remain empty and should be appropriately resized by the client code. */
    template<typename... Sizes> explicit ArrayTable(Sizes... sizes) { resize(sizes...); }

    /** This function resizes the table so that it holds the specified number of items in each of
        the N dimensions. If the last dimension in the list is nonzero, the arrays holding the data
        are resized correspondingly and all values are set to zero. If the last dimension is zero,
        the arrays are emptied and should be appropriately resized by the client code. In any case,
        any values that were previously in the table are lost. */
    template<typename... Sizes, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N, Sizes...>()>>
    void resize(Sizes... sizes)
    {
        _sizes = {{static_cast<size_t>(sizes)...}};
        _rows.resize(
            std::accumulate(_sizes.begin(), _sizes.end() - 1, static_cast<size_t>(1), std::multiplies<size_t>()));
        for (auto& row : _rows) row.resize(_sizes[N - 1]);
    }

    /** This function sets all values in the table to zero, without changing the number of items.
        */
    void setToZero()
    {
        for (auto& row : _rows) std::fill(begin(row), end(row), 0.);
    }

    // ================== Accessing sizes and values ==================

public:
    /** This function returns the number of items in the dimension indicated by the specified
        zero-based index. For the last dimension, the function returns the number of items as it
        has been specified in the constructor or in the most recent resize() invocation, even if
        the underlying arrays have since been resized by cient code. */
    size_t size(size_t dim) const { return _sizes[dim]; }

    /** This function returns the number of columns in a row, i.e. the number of items in the last
        dimension, defined as the size of the first row in the table, or zero if the table is
        empty. */
    size_t rowSize() const { return !_rows.empty() ? _rows[0].size() : 0; }

    /** This function returns an estimate for the total number of items in the table, obtained by
        multiplying the total number of rows in the table by the size of the first row. In other
        words, the estimate assumes that all rows have the same size. */
    size_t size() const { return _rows.size() * rowSize(); }

    /** This function returns a writable reference to the value at the specified N indices. There
        is no range checking. Out-of-range index values cause unpredictable behavior. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N, Indices...>()>>
    double& operator()(Indices... indices)
    {
        std::array<size_t, N> indexes = {{static_cast<size_t>(indices)...}};
        size_t rowIndex = indexes[0];
        for (size_t k = 1; k != N - 1; ++k) rowIndex = rowIndex * _sizes[k] + indexes[k];
        size_t columnIndex = indexes[N - 1];
        return _rows[rowIndex][columnIndex];
    }

    /** This function returns a copy of the value at the specified N indices. There is no range
        checking. Out-of-range index values cause unpredictable behavior. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N, Indices...>()>>
    double operator()(Indices... indices) const
    {
        std::array<size_t, N> indexes = {{static_cast<size_t>(indices)...}};
        size_t rowIndex = indexes[0];
        for (size_t k = 1; k != N - 1; ++k) rowIndex = rowIndex * _sizes[k] + indexes[k];
        size_t columnIndex = indexes[N - 1];
        return _rows[rowIndex][columnIndex];
    }

    /** This function returns a writable reference to the row at the specified N-1 indices. There
        is no range checking. Out-of-range index values cause unpredictable behavior. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N - 1, Indices...>()>>
    Array& operator()(Indices... indices)
    {
        std::array<size_t, N - 1> indexes = {{static_cast<size_t>(indices)...}};
        size_t rowIndex = indexes[0];
        for (size_t k = 1; k != N - 1; ++k) rowIndex = rowIndex * _sizes[k] + indexes[k];
        return _rows[rowIndex];
    }

    /** This function returns a read-only reference to the row at the specified N-1 indices. There
        is no range checking. Out-of-range index values cause unpredictable behavior. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N - 1, Indices...>()>>
    const Array& operator()(Indices... indices) const
    {
        std::array<size_t, N - 1> indexes = {{static_cast<size_t>(indices)...}};
        size_t rowIndex = indexes[0];
        for (size_t k = 1; k != N - 1; ++k) rowIndex = rowIndex * _sizes[k] + indexes[k];
        return _rows[rowIndex];
    }

    /** This function returns a writable reference to the row at the specified index (for a
        2-dimensional table only). There is no range checking. Out-of-range index values cause
        unpredictable behavior. */
    template<typename Index, typename = std::enable_if_t<N == 2 && CompileTimeUtils::isIntegralArgList<1, Index>()>>
    Array& operator[](Index index)
    {
        return _rows[index];
    }

    /** This function returns a read-only reference to the row at the specified index (for a
        2-dimensional table only). There is no range checking. Out-of-range index values cause
        unpredictable behavior. */
    template<typename Index, typename = std::enable_if_t<N == 2 && CompileTimeUtils::isIntegralArgList<1, Index>()>>
    const Array& operator[](Index index) const
    {
        return _rows[index];
    }

    // ================== Data members ==================

private:
    std::vector<Array> _rows;
    std::array<size_t, N> _sizes;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABLE_HPP
#define TABLE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include <array>
#include <functional>

////////////////////////////////////////////////////////////////////

/** An instance of the Table<N> class template holds an N-dimensional table of double values,
    offering indexed access for reading and writing individual values. All values are stored in a
    single data array. Values with adjacent rightmost indices are stored next to each other. */
template<size_t N> class Table
{
    static_assert(N >= 1, "Table dimension must be at least 1");

    // ================== Constructing ==================

public:
    /** The default constructor constructs an empty table. */
    Table() { _sizes.fill(0); }

    /** This constructor constructs a table holding the specified number of items in each of the N
        dimensions. All values are set to zero. */
    template<typename... Sizes> explicit Table(Sizes... sizes) { resize(sizes...); }

    /** This function resizes the table so that it holds the specified number of items in each of
        the N dimensions. All values are set to zero, i.e. any values that were previously in the
        table are lost. */
    template<typename... Sizes, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N, Sizes...>()>>
    void resize(Sizes... sizes)
    {
        _sizes = {{static_cast<size_t>(sizes)...}};
        _data.resize(std::accumulate(_sizes.begin(), _sizes.end(), static_cast<size_t>(1), std::multiplies<size_t>()));
    }

    /** This function sets all values in the table to zero, without changing the number of items.
        */
    void setToZero() { std::fill(begin(_data), end(_data), 0.); }

    // ================== Accessing sizes and values ==================

public:
    /** This function returns the total number of items in the table. */
    size_t size() const { return _data.size(); }

    /** This function returns the number of items in the dimension indicated by the specified
        zero-based index. */
    size_t size(size_t dim) const { return _sizes[dim]; }

    /** This function returns a writable reference to the value at the specified N indices. There
        is no range checking. Out-of-range index values cause unpredictable behavior. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N, Indices...>()>>
    double& operator()(Indices... indices)
    {
        return _data[flattenedIndex(indices...)];
    }

    /** This function returns a copy of the value at the specified N indices. There is no range
        checking. Out-of-range index values cause unpredictable behavior. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N, Indices...>()>>
    double operator()(Indices... indices) const
    {
        return _data[flattenedIndex(indices...)];
    }

    /** This function returns a writable reference to the value at the specified index (for a
        one-dimensional table only). There is no range checking. Out-of-range index values cause
        unpredictable behavior. */
    template<typename Index, typename = std::enable_if_t<N == 1 && CompileTimeUtils::isIntegralArgList<1, Index>()>>
    double& operator[](Index index)
    {
        return _data[index];
    }

    /** This function returns a copy of the value at the specified index (for a one-dimensional
        table only). There is no range checking. Out-of-range index values cause unpredictable
        behavior. */
    template<typename Index, typename = std::enable_if_t<N == 1 && CompileTimeUtils::isIntegralArgList<1, Index>()>>
    double operator[](Index index) const
    {
        return _data[index];
    }

    // ================== Accessing the raw data ==================

    /** This template function returns the flattened index in the underlying data array for the
        specified N indices. */
    template<typename... Indices, typename = std::enable_if_t<CompileTimeUtils::isIntegralArgList<N, Indices...>()>>
    size_t flattenedIndex(Indices... indices) const
    {
        std::array<size_t, N> indexes = {{static_cast<size_t>(indices)...}};
        size_t result = indexes[0];
        for (size_t k = 1; k != N; ++k) result = result * _sizes[k] + indexes[k];
        return result;
    }

    /** This function returns a writable reference to the underlying data array. Resizing the array
        other than through the Table::resize() function results in unpredictable behavior. */
    Array& data() { return _data; }

    /** This function returns a read-only reference to the underlying data array. */
    const Array& data() const { return _data; }

    // ================== Data members ==================

private:
    Array _data;
    std::array<size_t, N> _sizes;
};

////////////////////////////////////////////////////////////////////

#endif

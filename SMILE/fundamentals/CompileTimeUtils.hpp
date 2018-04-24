/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABLEIMPL_HPP
#define TABLEIMPL_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

// private compile-time utilities for Table and ArrayTable
namespace CompileTimeUtils
{
    // recursive implementation of "all" algorithm
    constexpr bool all() { return true; }
    template <typename... Ts>
    constexpr bool all(bool b, Ts... bs) { return b && all(bs...); }

    // returns true if the parameter pack has only integral arguments
    template <typename... Ts>
    inline constexpr bool isIntegralArgList()
    {
        return all(std::is_integral<Ts>::value...);
    }

    // returns true if the parameter pack has only floating point arguments
    template <typename... Ts>
    inline constexpr bool isFloatArgList()
    {
        return all(std::is_floating_point<Ts>::value...);
    }

    // returns true if the parameter pack has only integral arguments
    template <typename... Ts>
    inline constexpr bool isNumericArgList()
    {
        return all(std::is_arithmetic<Ts>::value...);
    }

    // returns true if the parameter pack has the specified number of integral arguments
    template <size_t N, typename... Ts>
    inline constexpr bool isIntegralArgList()
    {
        return sizeof...(Ts) == N && all(std::is_integral<Ts>::value...);
    }

    // returns true if the parameter pack has the specified number of floating point arguments
    template <size_t N, typename... Ts>
    inline constexpr bool isFloatArgList()
    {
        return sizeof...(Ts) == N && all(std::is_floating_point<Ts>::value...);
    }

    // returns true if the parameter pack has the specified number of integral arguments
    template <size_t N, typename... Ts>
    inline constexpr bool isNumericArgList()
    {
        return sizeof...(Ts) == N && all(std::is_arithmetic<Ts>::value...);
    }
}

////////////////////////////////////////////////////////////////////

#endif

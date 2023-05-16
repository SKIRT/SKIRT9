/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BASICS_HPP
#define BASICS_HPP

////////////////////////////////////////////////////////////////////

// Define this macro before including the system headers to silence warnings on Windows
// about C library functions that Microsoft considers to be unsafe, such as getenv()
#define _CRT_SECURE_NO_WARNINGS

// Define this macro before including the <cmath> system header to force the Windows
// version to define the de-facto standard mathematical constants such as M_PI
#define _USE_MATH_DEFINES

// Include the basic headers
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// Copy selected symbols to the global namespace
using std::abs;
using std::max;
using std::min;
using std::size_t;
using std::string;
using std::vector;

////////////////////////////////////////////////////////////////////

/** The Basics header file should be included directly or indirectly in all program units for the project
    to include the headers for some frequently-used standard library facilities in a consistent manner.

    The following standard library headers are included by this header:

    | header | description
    |--------|------------
    | \<algorithm\>        | general algorithms, e.g. all_of, count, find, sort
    | \<cfloat\>           | C-style floating point limits, e.g. DBL_MAX
    | \<climits\>          | C-style integral limits, e.g. INT_MAX
    | \<cmath\>            | mathematical functions, e.g. sin, exp, pow
    | \<cstddef\>          | standard types, e.g. size_t
    | \<cstdlib\>          | c-style functions, e.g. atol, strtol, div
    | \<initializer_list\> | initializer lists
    | \<iterator\>         | iterators and iterator support
    | \<limits\>           | C++ style limits, e.g. std::numeric_limits<double>::infinity()
    | \<memory\>           | resource management pointers
    | \<numeric\>          | numeric algorithms, e.g. accumulate
    | \<string\>           | sequence of characters
    | \<tuple\>            | tuples, tie
    | \<type_traits\>      | type traits
    | \<utility\>          | pair, swap, rvalue casts
    | \<vector\>           | one-dimensional resizable array

    In addition, the std::size_t, std::string, and std::vector types are copied to the global namespace
    so that the "std" prefix can be omitted for these types.

    Furthermore, the functions std::abs, std::min and std::max are copied to the global namespace
    so that the proper overloaded form is always used rather than the old C-style form. This is
    especially important for std::abs because the old C-style form converts doubles to integers.
*/
class Basics final
{};

////////////////////////////////////////////////////////////////////

#endif

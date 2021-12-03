/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NR_HPP
#define NR_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include "Range.hpp"

////////////////////////////////////////////////////////////////////

/** This static class contains a collection of functions that operate on Array objects and on
    std::vector<T> objects where T is some numeric type. Most implementations are provided inline
    in the header. */
class NR
{
    //======================== Conversion and Assignment =======================

public:
    /** This template function converts the source sequence to an Array object and returns the
        result. This template function works for any source container type with indexed random
        access (specifically including std::vector), and for any item types, as long as the source
        items can be assigned to or converted to double. */
    template<class V> static inline Array array(const V& sourcev)
    {
        size_t n = sourcev.size();
        Array resultv(n);
        for (size_t i = 0; i < n; i++) resultv[i] = sourcev[i];
        return resultv;
    }

    /** This template function assigns the source sequence to the destination array, resizing the
        destination array if necessary. This template function works for any source container type
        with indexed random access (specifically including std::vector), and for any item types, as
        long as the source items can be assigned to or converted to double. */
    template<class V, typename = std::enable_if_t<!CompileTimeUtils::isNumericArgList<V>()>>
    static inline void assign(Array& targetv, const V& sourcev)
    {
        size_t n = sourcev.size();
        targetv.resize(n);
        for (size_t i = 0; i < n; i++) targetv[i] = sourcev[i];
    }

    /** This template function assigns the floating point values specified as its arguments to the
        destination array, resizing the destination array if necessary. */
    template<typename... Values,
             typename = std::enable_if_t<sizeof...(Values) != 0 && CompileTimeUtils::isFloatArgList<Values...>()>>
    static inline void assign(Array& targetv, Values... values)
    {
        targetv.resize(sizeof...(Values));
        assignValues(0, targetv, values...);
    }

private:
    // recursively assign values from double arguments to Array
    static inline void assignValues(size_t /*index*/, Array& /*target*/) {}
    template<typename... Values>
    static inline void assignValues(size_t index, Array& target, double value, Values... values)
    {
        target[index] = value;
        assignValues(index + 1, target, values...);
    }

    //======================== Sorting =======================

public:
    /** This function sorts the values in the specified array. */
    static inline void sort(Array& xv) { std::sort(begin(xv), end(xv)); }

    /** This template function sorts a sequence of items given as a std::vector<T> where T is any
        built-in or user-defined type that implements the less-than operator (including the
        standard numeric types). */
    template<typename T> static inline void sort(std::vector<T>& xv) { std::sort(xv.begin(), xv.end()); }

    /** This template function sorts and removes duplicates from a sequence of items given as a
        std::vector<T> where T is any built-in or user-defined type that implements the less-than
        operator (including the standard numeric types). */
    template<typename T> static inline void unique(std::vector<T>& xv)
    {
        std::sort(xv.begin(), xv.end());
        size_t n = std::unique(xv.begin(), xv.end()) - xv.begin();
        xv.resize(n);
    }

    /** This function reverses blocks of values in the specified array. The number of values in the
        array must be an integer multiple of the specified block size. If not, the behavior is
        undefined. */
    static void reverse(Array& xv, size_t blocksize = 1);

    //======================== Searching =======================

public:
    /** This template function quickly performs a binary search on an ordered sequence of items
        provided as a std::vector<T>. It works for any item type T that implements the less-than
        operator (including the built-in numeric types). Given a sequence \f$\{x_i,\,i=0...N-1\}\f$
        and a query value \f$x\f$, the function returns the integer number \f$j\f$ such that \f$x_j
        \leq x < x_{j+1}\f$, as long as \f$x_0 \leq x < x_{N-1}\f$. In addition, if \f$x =
        x_{N-1}\f$ the function returns \f$N-2\f$; in other words the rightmost border is
        considered to be inside the last bin. If \f$x<x_0\f$ the function returns \f$-1\f$; if
        \f$x>x_{N-1}\f$ the function returns \f$N-1\f$. The function assumes that the provided
        sequence contains at least one element, and that the elements are sorted in ascending
        order. If this is not the case, the result is undefined. The algorithm is adapted from the
        Numerical Recipes in C++ handbook. */
    template<typename T> static inline int locate(const std::vector<T>& xv, const T& x)
    {
        int n = xv.size();
        if (x < xv[0]) return -1;
        if (xv[n - 1] < x) return n - 1;

        int jl = -1;
        int ju = n;
        while (ju - jl > 1)
        {
            int jm = (ju + jl) >> 1;
            if (x < xv[jm])
                ju = jm;
            else
                jl = jm;
        }
        if (jl <= 0) return 0;
        if (jl >= n - 2) return n - 2;
        return jl;
    }

    /** This function is an implementation detail and not intended for public use. See the
        description of the locate() function for more information. */
    static inline int locateBasicImpl(const Array& xv, double x, int n)
    {
        int jl = -1;
        int ju = n;
        while (ju - jl > 1)
        {
            int jm = (ju + jl) >> 1;
            if (x < xv[jm])
                ju = jm;
            else
                jl = jm;
        }
        return jl;
    }

    /** The locate(), locateClip() and locateFail() functions perform a binary search on the
        ordered sequence of double values in an array. There are subtle differences between the
        various "locate" functions, as decribed below. All functions assume that the specified
        array contains at least two elements, and that all elements are sorted in ascending order.
        If this is not the case, the result is undefined. The algorithm is adapted from the
        Numerical Recipes in C++ handbook. The array passed as the first function argument
        specifies an ordered sequence \f$\{x_i,\,i=0...N\}\f$ of \f$N+1\f$ values \f$x_i\f$,
        interpreted as borders defining a grid with \f$N\f$ bins. There must be at least two values
        (defining a single bin). The second function argument specifies a query value \f$x\f$. As
        long as the query is inside the range of the sequence, all functions return the zero-based
        index of the bin in which the query falls. Specifically, if \f$x_0 \leq x < x_N\f$ then all
        functions return the integer number \f$j\f$ such that \f$x_j \leq x < x_{j+1}\f$. In
        addition, if \f$x = x_N\f$ the function returns \f$N-1\f$. In other words the rightmost
        border is always considered to be inside the last bin (and of course the leftmost border is
        inside the first bin). The return value of the various "locate" functions differs for query
        values outside of the range of the sequence, as listed in the table below. <TABLE>
        <TR><TD><B>Function</B></TD> <TD><B>\f$x<x_0\f$</B></TD> <TD><B>\f$x>x_N\f$</B></TD>
        <TD><B>Comments</B> (note that there are \f$N+1\f$ values)</TD></TR> <TR><TD>locate()</TD>
        <TD>\f$-1\f$</TD> <TD>\f$N\f$</TD> <TD>Out-of-range values are indicated with a
        corresponding out-of-range index</TD></TR> <TR><TD>locateClip()</TD> <TD>\f$0\f$</TD>
        <TD>\f$N-1\f$</TD> <TD>Out-of-range values are considered to be inside the corresponding
        outermost bin</TD></TR> <TR><TD>locateFail()</TD> <TD>\f$-1\f$</TD> <TD>\f$-1\f$</TD>
        <TD>Out-of-range values are indicated with a negative index</TD></TR> </TABLE> */
    static inline int locate(const Array& xv, double x)
    {
        int n = xv.size();
        if (x == xv[n - 1]) return n - 2;
        return locateBasicImpl(xv, x, n);
    }

    /** This function performs a binary search on the ordered sequence of double values in an
        array. See the description of the locate() function for more information. */
    static inline int locateClip(const Array& xv, double x)
    {
        int n = xv.size();
        if (x < xv[0]) return 0;
        return locateBasicImpl(xv, x, n - 1);
    }

    /** This function performs a binary search on the ordered sequence of double values in an
        array. See the description of the locate() function for more information. */
    static inline int locateFail(const Array& xv, double x)
    {
        int n = xv.size();
        if (x > xv[n - 1]) return -1;
        return locateBasicImpl(xv, x, n - 1);
    }

    //======================== Constructing grids =======================

public:
    /** This function builds a linear grid over the specified range \f$[x_{\text{min}},
        x_{\text{max}}]\f$ and with the specified number of \f$N>0\f$ bins, and stores the
        resulting \f$N+1\f$ border points \f$x_i\f$ in the provided array, which is resized
        appropriately. The grid's border points are calculated according to \f[ x_i =
        x_{\text{min}} + (x_{\text{max}}-x_{\text{min}})\,\frac{i}{N} \qquad i=0,\ldots,N. \f] The
        function returns the bin width, i.e. the distance between two adjacent border points, given
        by \f$(x_{\text{max}}-x_{\text{min}})/N\f$. */
    static inline double buildLinearGrid(Array& xv, double xmin, double xmax, int n)
    {
        xv.resize(n + 1);
        double dx = (xmax - xmin) / n;
        for (int i = 0; i <= n; i++) xv[i] = xmin + i * dx;
        return dx;
    }

    /** This function builds a power-law grid over the specified range \f$[x_{\text{min}},
        x_{\text{max}}]\f$ and with the specified number of \f$N>0\f$ bins and specified ratio
        \f${\cal{R}}\f$ of last over first bin widths, and stores the resulting \f$N+1\f$ border
        points \f$x_i\f$ in the provided array, which is resized appropriately. If the specified
        ratio is very close to one, the function simply builds a linear grid. Otherwise the grid's
        border points are calculated according to \f[ x_i = x_{\text{min}} + (x_{\text{max}} -
        x_{\text{min}}) \left(\frac{1-q^i}{1-q^N}\right) \qquad i=0,\ldots,N, \f] with \f$ q =
        {\cal{R}}^{1/(N-1)} \f$. It is easy to check that the ratio between the widths of the last
        and first bins is indeed \f${\cal{R}}\f$: \f[ \frac{ x_N-x_{N-1} }{ x_1-x_0 } = \frac{
        q^{N-1} - q^N }{ 1-q } = q^{N-1} = {\cal{R}}. \f] */
    static inline void buildPowerLawGrid(Array& xv, double xmin, double xmax, int n, double ratio)
    {
        if (fabs(ratio - 1.) < 1e-3)
        {
            buildLinearGrid(xv, xmin, xmax, n);
        }
        else
        {
            xv.resize(n + 1);
            double range = xmax - xmin;
            double q = pow(ratio, 1. / (n - 1));
            double qn = pow(q, n);
            for (int i = 0; i <= n; ++i) xv[i] = xmin + (1. - pow(q, i)) / (1. - qn) * range;
        }
    }

    /** This function builds a symmetrical power-law grid over the specified range
        \f$[x_{\text{min}}, x_{\text{max}}]\f$ and with the specified number of \f$N>0\f$ bins and
        specified ratio \f${\cal{R}}\f$ of outermost over innermost bin widths, and stores the
        resulting \f$N+1\f$ border points \f$x_i\f$ in the provided array, which is resized
        appropriately. If the specified ratio is very close to one, the function simply builds a
        linear grid. Otherwise, because of the required symmetry, the expression for the grid's
        border points depends on whether the number of bins is even or odd. Define the centre of
        the grid as \f$x_{\text{c}} = \tfrac12(x_{\text{min}}+x_{\text{max}})\f$. If \f$N\f$ is
        even, we define \f$M=N/2\f$ and set \f$x_M=x_{\text{c}}\f$ and \f[ x_{M\pm i} =
        x_{\text{c}} \pm \tfrac12(x_{\text{max}}-x_{\text{min}}) \left(\frac{1-q^i}{1-q^M}\right)
        \qquad i=1,\ldots,M \f] with \f$ q = {\cal{R}}^{1/(M-1)} \f$. The ratio between the widths
        of the outermost and the innermost bins is now \f[ \frac{x_{2M}-x_{2M-1}}{x_{M+1}-x_M} =
        \frac{q^{M-1}-q^M}{1-q} = q^{M-1} = {\cal{R}}. \f] On the other hand, if \f$N\f$ is odd, we
        define \f$M = (N+1)/2\f$ and \f[ x_{M-\frac12\pm (i-\frac12)} = x_{\text{c}} \pm
        \tfrac12(x_{\text{max}}-x_{\text{min}}) \left[ \frac{ \frac12\,(1+q) - q^i }{
        \frac12\,(1+q) - q^M } \right] \qquad i=1,\ldots,M, \f] with again \f$ q =
        {\cal{R}}^{1/(M-1)} \f$. The ratio between the widths of the outermost and the innermost
        bins is for this case \f[ \frac{ x_{2M-1}-x_{2M-2} }{ x_M-x_{M-1} } = \frac{ q^{M-1} - q^M
        }{ 1-q } = q^{M-1} = {\cal{R}}. \f] */
    static inline void buildSymmetricPowerLawGrid(Array& xv, double xmin, double xmax, int n, double ratio)
    {
        if (fabs(ratio - 1.) < 1e-3)
        {
            buildLinearGrid(xv, xmin, xmax, n);
        }
        else
        {
            xv.resize(n + 1);
            double xc = 0.5 * (xmin + xmax);
            if (n % 2 == 0)
            {
                int M = n / 2;
                double q = pow(ratio, 1.0 / (M - 1.0));
                double qM = pow(q, M);
                xv[M] = xc;
                for (int i = 1; i <= M; ++i)
                {
                    double dxi = (1.0 - pow(q, i)) / (1.0 - qM) * 0.5 * (xmax - xmin);
                    xv[M + i] = xc + dxi;
                    xv[M - i] = xc - dxi;
                }
            }
            else
            {
                int M = (n + 1) / 2;
                double q = pow(ratio, 1.0 / (M - 1.0));
                double qM = pow(q, M);
                for (int i = 1; i <= M; ++i)
                {
                    double dxi = (0.5 + 0.5 * q - pow(q, i)) / (0.5 + 0.5 * q - qM) * 0.5 * (xmax - xmin);
                    xv[M - 1 + i] = xc + dxi;
                    xv[M - i] = xc - dxi;
                }
            }
        }
    }

    /** This function builds a logarithmic grid over the specified range \f$[x_{\text{min}},
        x_{\text{max}}]\f$ and with the specified number of \f$N>0\f$ bins, and stores the
        resulting \f$N+1\f$ border points \f$x_i\f$ in the provided array, which is resized
        appropriately. The grid's border points are calculated according to \f[ x_i =
        x_{\text{min}} \left( \frac{x_{\text{max}}}{x_{\text{min}}} \right)^{i/N} \qquad
        i=0,\ldots,N. \f] */
    static inline void buildLogGrid(Array& xv, double xmin, double xmax, int n)
    {
        xv.resize(n + 1);
        double logxmin = log(xmin);
        double dlogx = log(xmax / xmin) / n;
        for (int i = 0; i <= n; i++) xv[i] = exp(logxmin + i * dlogx);
    }

    /** This function builds a grid with its first bin starting at zero, and subsequent logarithmic
        border points over the specified range \f$[x_{\text{min}}, x_{\text{max}}]\f$. The grid has
        the specified number of \f$N>0\f$ bins and the resulting \f$N+1\f$ border points \f$x_i\f$
        are stored in the provided array, which is resized appropriately. The grid's border points
        are calculated according to \f$x_0=0\f$ and \f[ x_i = x_{\text{min}} \left(
        \frac{x_{\text{max}}}{x_{\text{min}}} \right)^{(i-1)/(N-1)} \qquad i=1,\ldots,N. \f] */
    static inline void buildZeroLogGrid(Array& xv, double xmin, double xmax, int n)
    {
        xv.resize(n + 1);  // this also initializes xv[0] to zero
        double logxmin = log(xmin);
        double dlogx = log(xmax / xmin) / (n - 1);
        for (int i = 0; i < n; i++) xv[i + 1] = exp(logxmin + i * dlogx);
    }

    //=================== Interpolating and resampling ===================

public:
    /** This function computes the interpolated value of a one-dimensional function, given its
        values at the edges of an interval. Both the axes coordinate \f$x\f$ and the function value
        \f$f(x)\f$ are interpolated linearly. */
    static inline double interpolateLinLin(double x, double x1, double x2, double f1, double f2)
    {
        return f1 + ((x - x1) / (x2 - x1)) * (f2 - f1);
    }

    /** This function computes the interpolated value of a one-dimensional function, given its
        values at the edges of an interval. The axes coordinate \f$x\f$ is interpolated
        logarithmically, and thus \f$x_1\f$ and \f$x_2\f$ must have positive values. If this is not
        the case, the function returns zero. The function value \f$f(x)\f$ is interpolated
        linearly. */
    static inline double interpolateLogLin(double x, double x1, double x2, double f1, double f2)
    {
        if (x1 <= 0 || x2 <= 0) return 0;
        return f1 + log(x / x1) / log(x2 / x1) * (f2 - f1);
    }

    /** This function computes the interpolated value of a one-dimensional function, given its
        values at the edges of an interval. The coordinates on both axes are interpolated
        logarithmically, and thus all arguments must have positive values. If the function values
        \f$f(x_1)\f$ and/or \f$f(x_2)\f$ are not positive, the function returns zero. The other
        values are not verified. */
    static inline double interpolateLogLog(double x, double x1, double x2, double f1, double f2)
    {
        if (f1 <= 0 || f2 <= 0)
        {
            if (x == x1) return f1;
            if (x == x2) return f2;
            return 0;
        }
        return f1 * exp(log(x / x1) / log(x2 / x1) * (log(f2 / f1)));
    }

    /** This template function returns the interpolated function value \f$y(x)\f$ given a value
        \f$x\f$ and a tabulated version of the function \f$y_k\f$ defined on a grid \f$x_k\f$ (both
        specified as arrays of the same length). If \f$x\f$ falls inside the grid, the \f$y(x)\f$
        value is interpolated between the neighboring grid points using the function specified as
        template argument. The interpolation functions provided by this namespace can be passed as
        a template argument for this purpose.

        If \f$x\f$ falls outside of the \f$x_k\f$ grid, zero is returned as the result.

        The function assumes that the specified arrays contain the same number of elements, that
        they each contain at least two elements, and that \f$x_k\f$ elements are sorted in
        ascending order. If this is not the case, the result is undefined. */
    template<double interpolateFunction(double, double, double, double, double)>
    static inline double value(double x, const Array& xv, const Array& yv)
    {
        int i = NR::locateFail(xv, x);
        if (i < 0) return 0.;
        return interpolateFunction(x, xv[i], xv[i + 1], yv[i], yv[i + 1]);
    }

    /** This template function returns the interpolated function value \f$y(x)\f$ given a value
        \f$x\f$ and a tabulated version of the function \f$y_k\f$ defined on a grid \f$x_k\f$ (both
        specified as arrays of the same length). If \f$x\f$ falls inside the grid, the \f$y(x)\f$
        value is interpolated between the neighboring grid points using the function specified as
        template argument. The interpolation functions provided by this namespace can be passed as
        a template argument for this purpose.

        If \f$x\f$ falls outside of the \f$x_k\f$ grid, the \f$y_k\f$ value at the corresponding
        outer grid point, i.e. \f$y_0\f$ or \f$y_{n-1}\f$ is returned as the result.

        The function assumes that the specified arrays contain the same number of elements, that
        they each contain at least two elements, and that \f$x_k\f$ elements are sorted in
        ascending order. If this is not the case, the result is undefined. */
    template<double interpolateFunction(double, double, double, double, double)>
    static inline double clampedValue(double x, const Array& xv, const Array& yv)
    {
        int n = xv.size();
        int i = NR::locate(xv, x);
        if (i < 0) return yv[0];
        if (i >= n - 1) return yv[n - 1];
        return interpolateFunction(x, xv[i], xv[i + 1], yv[i], yv[i + 1]);
    }

    /** This template function resamples the function values \f$y_k\f$ defined on a grid \f$x_k\f$
        (both specified as arrays of the same length) onto an new grid \f$x_l^*\f$. The result is
        returned as an array of function values \f$y_l^*\f$ with the same length as the target
        grid. For new grid points that fall beyond the original grid, the function value is set to
        zero. For new grid points inside the original grid, the function value is interpolated
        using the function specified as template argument. The interpolation functions provided by
        this namespace can be passed as a template argument for this purpose. */
    template<double interpolateFunction(double, double, double, double, double)>
    static inline Array resample(const Array& xresv, const Array& xoriv, const Array& yoriv)
    {
        int n = xresv.size();
        Array yresv(n);
        for (int l = 0; l < n; l++) yresv[l] = value<interpolateFunction>(xresv[l], xoriv, yoriv);
        return yresv;
    }

    /** This template function resamples the function values \f$y_k\f$ defined on a grid \f$x_k\f$
        (both specified as arrays of the same length) onto an new grid \f$x_l^*\f$. The result is
        returned as an array of function values \f$y_l^*\f$ with the same length as the target
        grid. For new grid points that fall beyond the original grid, the function value is set to
        the value at the corresponding outer grid point, i.e. \f$y_0\f$ or \f$y_{n-1}\f$. For new
        grid points inside the original grid, the function value is interpolated using the function
        specified as template argument. The interpolation functions provided by this namespace can
        be passed as a template argument for this purpose. */
    template<double interpolateFunction(double, double, double, double, double)>
    static inline Array clampedResample(const Array& xresv, const Array& xoriv, const Array& yoriv)
    {
        int n = xresv.size();
        Array yresv(n);
        for (int l = 0; l < n; l++) yresv[l] = clampedValue<interpolateFunction>(xresv[l], xoriv, yoriv);
        return yresv;
    }

    //=============== Constructing cumulative distribution functions ==================

public:
    /** Given a discrete distribution over \f$N\f$ points \f[p_i \qquad i=0,\dots,N-1\f] this
        function builds the corresponding normalized cumulative distribution \f[P_0=0;\quad
        P_{i+1}=\frac{\sum_{j=0}^i p_j}{\sum_{j=0}^{N-1} p_j} \qquad i=0,\dots,N-1\f] with
        \f$N+1\f$ elements. The function returns the factor used to normalize the cumulative
        distribution. In this version of the function, the source distribution is specified as an
        array with at least one element; the target array is resized appropriately and replaced by
        the cumulative distribution. */
    static inline double cdf(Array& Pv, const Array& pv)
    {
        int n = pv.size();
        Pv.resize(n + 1);  // also sets Pv[0] to zero
        for (int i = 0; i < n; i++) Pv[i + 1] = Pv[i] + pv[i];
        double norm = Pv[n];
        Pv /= norm;
        return norm;
    }

    /** Given a discrete distribution over \f$N\f$ points \f[p_i \qquad i=0,\dots,N-1\f] this
        function builds the corresponding normalized cumulative distribution \f[P_0=0;\quad
        P_{i+1}=\frac{\sum_{j=0}^i p_j}{\sum_{j=0}^{N-1} p_j} \qquad i=0,\dots,N-1\f] with
        \f$N+1\f$ elements. The function returns the factor used to normalize the cumulative
        distribution. In this version of the function, the source distribution is specified by a
        function object with signature double pv(int i); the number of source points \f$N>0\f$ is
        specified as a separate argument. The source function is called once for each index
        \f$i=0,\dots,N-1\f$. The target array is resized appropriately and replaced by the
        cumulative distribution. */
    template<typename Functor> static inline double cdf(Array& Pv, int n, Functor pv)
    {
        Pv.resize(n + 1);  // also sets Pv[0] to zero
        for (int i = 0; i < n; i++) Pv[i + 1] = Pv[i] + pv(i);
        double norm = Pv[n];
        Pv /= norm;
        return norm;
    }

    /** Given the tabulated values for a continuous probability distribution, this function
        constructs both the normalized probability density function (pdf) and the corresponding
        normalized cumulative distribution function (cdf) within a given range. The incoming
        distribution is specified by the \em inxv and inpv arrays and does not need to be
        normalized. The specified \em xrange must overlap the incoming grid, but it does not need
        to coincide with it. The function constructs a new grid \em xv that matches the given range
        and then constructs both the normalized pdf \em pv and the normalized cdf \em Pv
        corresponding to this new grid. It returns the factor used to normalize the distributions.

        The outer grid points of the new grid match the specified range, and therefore usually do
        not occur in the incoming grid. The \em pv values for those new grid points are
        interpolated using the function specified as template argument. The interpolation functions
        provided by this namespace can be passed as a template argument for this purpose.

        The interpolation function specified as a template argument also determines the integration
        mechanism used to calculate the cdf. If the NR::interpolateLogLog() function is specified,
        it is assumed that the pdf behaves as a power-law between any two grid points, and the
        integration is performed accordingly. In all other cases, piece-wise linear behavior is
        assumed and regular trapezium-rule integration is used. */
    template<double interpolateFunction(double, double, double, double, double)>
    static inline double cdf(Array& xv, Array& pv, Array& Pv, const Array& inxv, const Array& inpv, Range xrange)
    {
        // copy the relevant portion of the axis grid
        size_t minRight = std::upper_bound(begin(inxv), end(inxv), xrange.min()) - begin(inxv);
        size_t maxRight = std::lower_bound(begin(inxv), end(inxv), xrange.max()) - begin(inxv);
        size_t n = 1 + maxRight - minRight;  // n = number of bins
        xv.resize(n + 1);                    // n+1 = number of border points
        size_t i = 0;                        // i = index in target array
        xv[i++] = xrange.min();              // j = index in input array
        for (size_t j = minRight; j < maxRight;) xv[i++] = inxv[j++];
        xv[i++] = xrange.max();

        // interpolate or copy the corresponding probability density values
        pv.resize(n + 1);
        pv[0] = minRight == 0 ? 0.
                              : interpolateFunction(xv[0], inxv[minRight - 1], inxv[minRight], inpv[minRight - 1],
                                                    inpv[minRight]);
        for (size_t i = 1; i < n; ++i) pv[i] = inpv[minRight + i - 1];
        pv[n] = maxRight == inxv.size() ? 0.
                                        : interpolateFunction(xv[n], inxv[maxRight - 1], inxv[maxRight],
                                                              inpv[maxRight - 1], inpv[maxRight]);

        // perform the rest of the operation in a non-templated function
        // assume log-log behavior for the cdf integration if the template parameter is NR::interpolateLogLog
        //  --> this will break if the user specifies another log-log interpolation function
        return cdf2(interpolateFunction == NR::interpolateLogLog, xv, pv, Pv);
    }

    /** This function calculates the normalized cdf \em Pv for the given axis grid \em xv and
        corresponding unnormalized pdf \em pv. It also normalizes the incoming pdf \em pv and
        returns the normalization factor. This non-template function is not implemented inline (1)
        to avoid code-duplication and (2) to avoid propagating the dependencies (includes) of the
        implementation to all NR users.

        If the \em loglog flag is false, piece-wise linear behavior of both the pdf and cdf is
        assumed and regular trapezium-rule integration is used. If the \em loglog flag is true, it
        is assumed that the pdf is linear in log-log space between any two grid points (equivalent
        to power-law behavior), and the integration is performed accordingly, as described below.

        Consider the pdf values \f$p_i\f$ and \f$p_{i+1}\f$ at two consecutive grid points
        \f$x_i\f$ and \f$x_{i+1}\f$. Assuming power-law behavior, the pdf between these two grid
        points can be written as \f[ p(x) = p_i \left(\frac{x}{x_i}\right)^{\alpha_i},
        \quad\mathrm{with}\; \alpha_i = \frac{\ln p_{i+1}/\ln p_i}{\ln x_{i+1}/\ln x_i} \f]

        The area under the curve is then \f[ \int_{x_i}^{x_{i+1}} p(x)\,\mathrm{d}x =
        \int_{x_i}^{x_{i+1}} p_i \left(\frac{x}{x_i}\right)^{\alpha_i}\mathrm{d}x = p_i x_i
        \;\mathrm{gln}\left(-\alpha_i, \frac{x_{i+1}}{x_i}\right) \f] where \f$\mathrm{gln}(a,x)\f$
        is the generalized logarithm defined in the description of the SpecialFunctions::gln()
        function. */
    static double cdf2(bool loglog, const Array& xv, Array& pv, Array& Pv);
};

////////////////////////////////////////////////////////////////////

#endif

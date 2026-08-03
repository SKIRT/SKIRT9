/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DISJOINTWAVELENGTHGRID_HPP
#define DISJOINTWAVELENGTHGRID_HPP

#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** DisjointWavelengthGrid is an abstract class that represents wavelength grids with
    straightforward, non-overlapping bins with constant transmission across each bin.

    Specifically, a disjoint wavelength grid consists of non-overlapping but possibly adjacent
    wavelength bins in increasing wavelength order, with constant maximum transmission within the
    bins and zero transmission outside of the bins. Each bin is defined by its left and right
    borders and has a characteristic wavelength that falls inside the bin. The left border is
    considered to be inside of the bin; the right border is considered to be outside of the bin.
    Neighboring bins may have a common border but can also be disconnected.

    Formally, assuming \f$N>0\f$ bins with zero-based indices, we have \f[
    \lambda^\mathrm{left}_\ell \le \lambda^\mathrm{c}_\ell < \lambda^\mathrm{right}_\ell, \quad
    \ell=0\dots N-1 \f] and if \f$N>1\f$, we additionally have \f[ \lambda^\mathrm{right}_\ell \le
    \lambda^\mathrm{left}_{\ell+1}, \quad \ell=0\dots N-2. \f] Finally, each bin of course has an
    associated bin width, \f[\lambda^\mathrm{right}_\ell - \lambda^\mathrm{left}_\ell > 0, \quad
    \ell=0\dots N-1.\f]

    A DisjointWavelengthGrid subclass is expected to invoke one of the setWavelengthXXX() functions
    during setup to initialize the wavelength grid. The current implementation offers two such
    functions: one to specify a consecutive range of adjacent wavelength bins given a list of
    characteric wavelengths, and another one to specify distinct, nonadjacent wavelength bins given
    a list of characteric wavelengths and a relative bin width. Other options can be added as the
    need arises. */
class DisjointWavelengthGrid : public WavelengthGrid
{
    ITEM_ABSTRACT(DisjointWavelengthGrid, WavelengthGrid, "a wavelength grid with non-overlapping bins")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the wavelength bins have been initialized by a subclass calling
        one of the setWavelengthXXX() functions of this class in their setupSelfBefore() function.
        */
    void setupSelfAfter() override;

    /** This function initializes the wavelength grid to a consecutive range of adjacent wavelength
        bins given a list of characteric wavelengths. The subclass determines a list of characteric
        wavelengths according to some predefined scheme, and the bin borders and bin widths are
        automatically determined from that list by this function. If the specified wavelength list
        is empty, or if there are duplicate values (which would lead to empty bins), the function
        throws a fatal error.

        The function first sorts the specified characteristic wavelengths in ascending order and
        then calculates the bin borders assuming linear or logarithmic scaling depending on the
        value of the \em logScale flag. The inner border between two consecutive characteristic
        wavelengths is placed at the mid-point of those two wavelengths (in linear or logarithmic
        space), and the outer borders at the edges of the wavelength range are placed such that the
        outer bins have the same width as the respective adjacent bins (in linear or logarithmic
        space). If the wavelength grid has just a single wavelength, the borders are placed just
        next to the wavelength to form a narrow bin. Finally, the function trivially calculates the
        wavelength bin widths from the bin borders.

        For linear scaling, the corresponding formulas are trivial. For logarithmic scaling, the
        formalas in logarithmic space translate easily to equivalent but more efficient formulas in
        real space. For the inner borders this yields the geometric mean of the two adjacent
        characteristic wavelengths, i.e. \f$\lambda^\mathrm{right}_{\ell-1} =
        \lambda^\mathrm{left}_\ell = \sqrt{\lambda^\mathrm{c}_{\ell-1}\lambda^\mathrm{c}_\ell}\;,
        \ell=1\dots N-1\f$. The leftmost outer border is placed at \f$\lambda^\mathrm{left}_0 =
        \sqrt{(\lambda^\mathrm{c}_{0})^3/\lambda^\mathrm{c}_1}\f$, and the rightmost outer border
        is placed at \f$\lambda^\mathrm{right}_{N-1} =
        \sqrt{(\lambda^\mathrm{c}_{N-1})^3/\lambda^\mathrm{c}_{N-2}}\f$.

        If there is just a single wavelength in the grid, the outer borders are placed (for both
        linear and logarithmic scaling) according to \f$\lambda^\mathrm{left}_0 =
        \lambda^\mathrm{c}_{0}(1-1/1000)\f$ and \f$\lambda^\mathrm{right}_0 =
        \lambda^\mathrm{c}_{0}(1+1/1000)\f$. */
    void setWavelengthRange(const Array& lambdav, bool logScale);

    /** This function initializes the wavelength grid to a set of distinct, nonadjacent wavelength
        bins given a list of characteric wavelengths and a relative half bin width. The subclass
        determines a list of characteric wavelengths and a relative half bin width, and the bin
        borders and bin widths are automatically calculated from that information by this function.
        If the specified wavelength list is empty, or if the relative half bin width is not
        positive, or if the calculated bins overlap, the function throws a fatal error.

        Specifically, the function first sorts the specified characteristic wavelengths in
        ascending order and then calculates the bin borders using \f$\lambda^\mathrm{left}_\ell =
        \lambda^\mathrm{c}_\ell(1-w)\f$ and \f$\lambda^\mathrm{right}_\ell =
        \lambda^\mathrm{c}_\ell(1+w)\;, \ell=0\dots N-1\f$, where \f$w\f$ is the specified relative
        half bin width. If the \em constantWidth flag is true, the width for the shortest
        wavelength is used for all bin widths instead. Finally the function trivially calculates
        the wavelength bin widths from the bin borders. */
    void setWavelengthBins(const Array& lambdav, double relativeHalfWidth, bool constantWidth = false);

    /** This function initializes the wavelength grid to a consecutive range of \f$N>0\f$ adjacent
        wavelength bins given a list of \f$N+1\f$ wavelength bin borders. The subclass determines a
        list of bin borders according to some predefined scheme, and the characteristic wavelengths
        and bin widths are automatically determined from that list by this function. If the
        specified list has fewer than two bin borders, or if there are duplicate values (which
        would lead to empty bins), the function throws a fatal error.

        The function first sorts the specified wavelength bin borders in ascending order and then
        calculates the characteristic wavelengths assuming linear scaling (arithmetic mean) or
        logarithmic scaling (geometric mean) depending on the value of the \em logScale flag. */
    void setWavelengthBorders(const Array& borderv, bool logScale);

    /** This function initializes the wavelength grid from a list of interleaved bin border points
        and corresponding characteristic wavelengths. (i.e., borders and characteristic wavelengths
        alternate). The number of values must be uneven and at least three. The list must be in
        strictly increasing or decreasing order, which means duplicates are not allowed, except
        that a zero characteristic wavelength indicates a segment that is not part of the grid,
        i.e. that lies between two non-adjacent bins. In other words, this option allows to (1)
        arbitrarily place characteristic wavelengths within each bin and (2) to specify
        intermediate wavelength ranges that are not covered by any bin. */
    void setWavelengthSegments(const Array& bordcharv);

    /** This function initializes the wavelength grid from a list of bin border points and a
        corresponding list of characteristic wavelengths. A zero characteristic wavelength
        indicates a segment that is not part of the grid, i.e. that lies between two non-adjacent
        bins or beyond the outer grid borders. In other words, the subclass has full control over
        the placement of bin borders and characteristic wavelengths.

        The two lists must have the same size. The borders must be listed in strictly increasing
        order of wavelength (i.e. there cannot be any duplicates), all nonzero characteristic
        wavelengths must lie within their corresponding bin, the last characteristic wavelength
        must be zero, and there must be at least one nonzero characteristic wavelength. If these
        requirements are violated, the behavior of this function is undefined. */
    void setWavelengthSegments(const vector<double>& borderv, const vector<double>& characv);

    //================= Functions implementing virtual base class functions ===================

public:
    /** This function returns the number of bins, \f$N\f$, in the grid (or equivalently, the number
        of characteristic wavelengths). */
    int numBins() const override;

    /** This function returns the characteristic wavelength \f$\lambda^\mathrm{c}_\ell\f$
       corresponding to the index \f$\ell\f$. */
    double wavelength(int ell) const override;

    /** This function returns the left border of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{left}_\ell\f$. */
    double leftBorder(int ell) const override;

    /** This function returns the right border of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{right}_\ell\f$. */
    double rightBorder(int ell) const override;

    /** This function returns the width of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{right}_\ell - \lambda^\mathrm{left}_\ell\f$. */
    double effectiveWidth(int ell) const override;

    /** This function returns the relative transmission for the wavelength bin corresponding to the
        index \f$\ell\f$ at the wavelength \f$\lambda\f$. For the present class, it always returns
        1. */
    double transmission(int ell, double lambda) const override;

    /** This function returns a single-element list with the index \f$\ell\f$ of the wavelength bin
        that contains the specified wavelength \f$\lambda\f$, i.e. for which
        \f$\lambda^\mathrm{left}_\ell <= \lambda < \lambda^\mathrm{right}_\ell\f$. If \f$\lambda\f$
        does not lie inside one of the wavelength bins, the function returns an empty list. */
    vector<int> bins(double lambda) const override;

    /** This function returns the index \f$\ell\f$ of the wavelength bin that contains the
        specified wavelength \f$\lambda\f$, i.e. for which \f$\lambda^\mathrm{left}_\ell <= \lambda
        < \lambda^\mathrm{right}_\ell\f$. If \f$\lambda\f$ does not lie inside one of the
        wavelength bins, the function returns -1. */
    int bin(double lambda) const override;

    //=============== Functions specific to disjoint wavelength grids =================

public:
    /** This function returns (a reference to) the list of characteristic wavelengths in this
        wavelength grid. In combination with the dlambdav() function, it allows easily expressing
        calculations involving consecutive wavelength grids. */
    const Array& lambdav() const { return _lambdav; }

    /** This function returns (a reference to) the list of bin widths in this wavelength grid. In
        combination with the lambdav() function, it allows easily expressing calculations involving
        consecutive wavelength grids. */
    const Array& dlambdav() const { return _dlambdav; }

    /** This function returns a list of the characteristic wavelengths in this wavelength grid
        extended with the outermost bin border point on each side. The list has thus two additional
        points, one on each side, and as a result covers the complete wavelength range of the grid,
        including the widths of the outer bins. This extended wavelength list can be used in
        situations where one needs to calculate/interpolate some function over the complete range
        of the wavelength grid and not just up to the outermost characteristic wavelengths. */
    Array extlambdav() const;

    /** This function returns a list of bin widths in this wavelength grid extended with a zero
        value on each side. This extended wavelength bin width list can be used, for example, to
        integrate a function discretized on the extended wavelength grid returned by the
        extlambdav() function over the wavelength range. */
    Array extdlambdav() const;

    //======================== Data Members ========================

private:
    // subclasses should call setWavelengthXXX() in their setupSelfBefore() to initialize these tables
    Array _lambdav;       // N characteristic wavelengths
    Array _dlambdav;      // N wavelength bin widths
    Array _lambdaleftv;   // N left wavelength bin borders
    Array _lambdarightv;  // N right wavelength bin widths
    Array _borderv;       // K=N+1 or K=N*2 ordered border points (depending on whether bins are adjacent)
    vector<int> _ellv;    // K+1 indices of the wavelength bins defined by the border points, or -1 if out of range
};

//////////////////////////////////////////////////////////////////////

#endif

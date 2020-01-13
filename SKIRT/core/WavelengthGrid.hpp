/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHGRID_HPP
#define WAVELENGTHGRID_HPP

#include "Array.hpp"
#include "Range.hpp"
#include "SimulationItem.hpp"

//////////////////////////////////////////////////////////////////////

/** WavelengthGrid is an abstract class that defines the interface for wavelength grids that can be
    used, for example, to specify the wavelength bins used for detetecting photon packets in
    instruments.

    A wavelength grid consists of \f$N>0\f$ possibly overlapping wavelength bins. Generally
    speaking, each of these bins is defined through a transmission curve determining the
    contribution fraction of a detected photon packet to the bin as a function of wavelength. This
    class defines the public interface for a wavelength grid in these general terms. Subclasses can
    implement several bin types, including bins that mimic a particular broadband, or
    straightforward bins with constant transmission across some wavelength interval.

    Key properties of a wavelength bin include its left and right borders, defining a wavelength
    interval within which the transmission may be nonzero, and its characteristic wavelength, which
    will be used by instruments to convert mean specific luminosities in the band between
    wavelength and frequency representations. The WavelengthGrid class requires that the
    characteristic wavelength of a bin falls inside the bin, i.e. \f$\lambda^\mathrm{left}_\ell \le
    \lambda^\mathrm{c}_\ell \le \lambda^\mathrm{right}_\ell, \ell=0\dots N-1\f$, and that no two bins
    have the same characteristic wavelength. This allows the WavelengthGrid class to meaningfully
    sort bins in increasing order of characteristic wavelength, and assign bin indices accordingly.

    The public interface also includes functions to obtain the relative transmission for a given
    bin as a function of wavelength, defined as the transmission at that wavelength divided by the
    maximum transmission, and to obtain a bin's effective width, defined as the horizontal size of
    a rectangle with height equal to the maximum transmission and with the same area as the one
    covered by the band's transmission curve.

    Finally, and most importantly, the public interface offers a function to determine the (indices
    of) the bin(s) that may have a nonzero transmission at a given wavelength. */
class WavelengthGrid : public SimulationItem
{
    ITEM_ABSTRACT(WavelengthGrid, SimulationItem, "a wavelength grid")
    ITEM_END()

    //======================== Public interface =======================

public:
    /** This function returns the number of bins \f$N\f$ in the grid, or equivalently, the number
        of characteristic wavelengths. Bins are always sorted in order of increasing characteristic
        wavelength. */
    virtual int numBins() const = 0;

    /** This function returns the characteristic wavelength \f$\lambda^\mathrm{c}_\ell\f$
       corresponding to the index \f$\ell\f$. The characteristic wavelength of a bin is always
       inside the bin, i.e. \f$\lambda^\mathrm{left}_\ell \le \lambda^\mathrm{c}_\ell \le
       \lambda^\mathrm{right}_\ell\f$. It can be used to convert mean specific luminosities in the
       band between wavelength and frequency representations. Bins are always sorted in order of
       increasing characteristic wavelength. */
    virtual double wavelength(int ell) const = 0;

    /** This function returns the left border of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{left}_\ell\f$. The transmission for this bin is
        guaranteed to be zero for all wavelengths shorter than the left border. */
    virtual double leftBorder(int ell) const = 0;

    /** This function returns the right border of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{right}_\ell\f$. The transmission for this bin is
        guaranteed to be zero for all wavelengths longer than the right border. */
    virtual double rightBorder(int ell) const = 0;

    /** This function returns the effective width of the wavelength bin corresponding to the index
        \f$\ell\f$. The effective width is defined as the horizontal size of a rectangle with
        height equal to the maximum transmission and with the same area as the one covered by the
        band's transmission curve. For bins with a constant transmission over the complete
        interval, the effective width is simply \f$\lambda^\mathrm{right}_\ell -
        \lambda^\mathrm{left}_\ell\f$. For bins with non-constant transmission, this is no longer
        true. */
    virtual double effectiveWidth(int ell) const = 0;

    /** This function returns the relative transmission for the wavelength bin corresponding to the
        index \f$\ell\f$ at the wavelength \f$\lambda\f$. The relative transmission is defined as
        the transmission at that wavelength divided by the maximum transmission for the bin. */
    virtual double transmission(int ell, double lambda) const = 0;

    /** This function returns a list of indices \f$\ell_k\f$ of the wavelength bins that may have a
        nonzero transmission at the specified wavelength \f$\lambda\f$, i.e. for which
        \f$\lambda^\mathrm{left}_\ell \le \lambda \le \lambda^\mathrm{right}_\ell\f$. If no
        wavelengths bins match this condition, the function returns an empty list. */
    virtual vector<int> bins(double lambda) const = 0;

    /** This function returns the index \f$\ell\f$ of one the wavelength bins that may have a
        nonzero transmission at the specified wavelength \f$\lambda\f$, i.e. for which
        \f$\lambda^\mathrm{left}_\ell \le \lambda \le \lambda^\mathrm{right}_\ell\f$. If no
        wavelengths bins match this condition, the function returns -1. If multiple bins match this
        condition, the function returns the index for the bin with the shortest characteristic
        wavelength. */
    virtual int bin(double lambda) const = 0;

    /** This function returns the wavelength range covered by the wavelength grid, which is defined
        as the range from the left border of the leftmost bin to the right border of the rightmost
        bin. This range includes all wavelengths possibly covered by the wavelength grid except in
        the rare case of overlapping bins where an inner bin is so wide that it outer limit extends
        beyond the outer bin. */
    Range wavelengthRange() const;
};

//////////////////////////////////////////////////////////////////////

#endif

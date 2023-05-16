/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BANDWAVELENGTHGRID_HPP
#define BANDWAVELENGTHGRID_HPP

#include "Band.hpp"
#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** BandWavelengthGrid is an asbtract class for representing wavelength grids where each bin is
    defined by a Band object. Each Band object provides the transmission curve and the other
    (derived) properties for the corresponding bin. Examples include the standard Johnson filters
    and the broadbands for actual observatories such as GALEX, SDSS or Herschel. See the Band class
    for more information.

    Each subclass must override the bandList() method to return a list of Band objects. The bands
    are automatically sorted in order of increasing pivot wavelength, even if the list returned by
    the subclass has them in a different order. The intervals in which the transmission is nonzero
    is allowed to overlap, but no two bands in the list should have the same pivot wavelength (this
    restriction essentially disallows specifying the same band twice). */
class BandWavelengthGrid : public WavelengthGrid
{
    ITEM_ABSTRACT(BandWavelengthGrid, WavelengthGrid, "a wavelength grid defined by a list of (broad)bands")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains the bands from the subclass and sorts them in order of increasing
        pivot wavelength. */
    void setupSelfAfter() override;

    /** This function must be implemented by each subclass to return a list of Band objects in
        arbitrary order (the bands will be sorted in order of increasing pivot wavelength by this
        class). */
    virtual vector<Band*> bandList() = 0;

    //================= Functions implementing virtual base class functions ===================

public:
    /** This function returns the number of bands in the grid (or equivalently, the number of
        bins). Bands are always sorted in order of increasing pivot wavelength. */
    int numBins() const override;

    /** This function returns the pivot wavelength for the band corresponding to the index
       \f$\ell\f$. Refer to the Band class for the relevant formulas. Bands are always sorted in
       order of increasing pivot wavelength. */
    double wavelength(int ell) const override;

    /** This function returns the left wavelength border of the band corresponding to the index
        \f$\ell\f$. The transmission for this band is guaranteed to be zero for all wavelengths
        shorter than the left border. */
    double leftBorder(int ell) const override;

    /** This function returns the right wavelength border of the band corresponding to the index
        \f$\ell\f$. The transmission for this band is guaranteed to be zero for all wavelengths
        longer than the right border. */
    double rightBorder(int ell) const override;

    /** This function returns the effective width of the band corresponding to the index
        \f$\ell\f$. Refer to the Band class for the relevant formulas. */
    double effectiveWidth(int ell) const override;

    /** This function returns the relative transmission for the band corresponding to the index
        \f$\ell\f$ at the wavelength \f$\lambda\f$. Refer to the Band class for the relevant
        formulas. */
    double transmission(int ell, double lambda) const override;

    /** This function returns a list of indices \f$\ell_k\f$ of the bands that may have a nonzero
        transmission at the specified wavelength \f$\lambda\f$, i.e. for which
        \f$\lambda^\mathrm{left}_\ell \le \lambda \le \lambda^\mathrm{right}_\ell\f$. If no bands
        match this condition, the function returns an empty list. */
    vector<int> bins(double lambda) const override;

    /** This function returns the index \f$\ell\f$ of a band that may have a nonzero transmission
        at the specified wavelength \f$\lambda\f$, i.e. for which \f$\lambda^\mathrm{left}_\ell \le
        \lambda \le \lambda^\mathrm{right}_\ell\f$. If no bands match this condition, the function
        returns -1. If multiple bands match this condition, the function returns the index for the
        band with the shortest characteristic wavelength. */
    int bin(double lambda) const override;

    //=============== Functions specific to band wavelength grids =================

public:
    /** This function returns (a pointer to) the Band object defining the transmission curve for
        the bin with index \f$\ell\f$. */
    const Band* band(int ell) const;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    vector<Band*> _bands;
};

//////////////////////////////////////////////////////////////////////

#endif

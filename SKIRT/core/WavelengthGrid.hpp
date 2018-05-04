/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHGRID_HPP
#define WAVELENGTHGRID_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"

//////////////////////////////////////////////////////////////////////

/** WavelengthGrid is an abstract class that defines the interface for wavelength grids that can be
    used, for example, to specify the wavelength discretization in instruments. A wavelength grid
    consists of \f$N>0\f$ consecutive wavelength bins \f$[\lambda^\mathrm{b}_\ell,
    \lambda^\mathrm{b}_{\ell+1}] \;,\ell=0\dots N-1\f$, defined by \f$N+1\f$ border points
    \f$\lambda^\mathrm{b}_\ell \;,\ell=0\dots N\f$. Each bin has an associated representative
    wavelength \f$\lambda_\ell \;,\ell=0\dots N-1\f$ that falls inside the bin without touching the
    borders, i.e. \f$\lambda^\mathrm{b}_\ell < \lambda_\ell < \lambda^\mathrm{b}_{\ell+1}
    \;,\ell=0\dots N-1\f$. Also, each bin of course has an associated bin width,
    \f$\lambda^\mathrm{b}_{\ell+1} - \lambda^\mathrm{b}_\ell \;,\ell=0\dots N-1\f$.

    In practice, one of the WavelengthGrid subclasses specifies the list of representative
    wavelengths according to some predefined scheme, and the bin borders and bin widths are
    automatically determined from this list. */
class WavelengthGrid : public SimulationItem
{
    ITEM_ABSTRACT(WavelengthGrid, SimulationItem, "a wavelength grid")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the wavelength bins have been initialized by a subclass calling
        the setWavelengths() function of this class in their setupSelfBefore() function. */
    void setupSelfAfter() override;

    /** This function sets the wavelength grid to the given list of representative wavelengths, and
        calculates the corresponding bin borders and bin widths. It should be called from the
        setupSelfBefore() function in each WavelengthGrid subclass. If the specified wavelength
        array is empty, of if there are duplicate values (which would lead to empty bins), the
        function throws a fatal error.

        The function first sorts the specified representative wavelengths in ascending order and
        then calculates the bin borders as follows. The inner border between two consecutive
        representative wavelengths is placed at the geometric mean of those two wavelengths, i.e.
        \f$\lambda^\mathrm{b}_\ell = \sqrt{\lambda_{\ell-1}\lambda_\ell} \;,\ell=1\dots N-1\f$.
        The outer borders at the very left and right of the wavelength range are placed just
        outside of the range, i.e. \f$\lambda^\mathrm{b}_0 = \lambda_0(1-1/1000)\f$ and
        \f$\lambda^\mathrm{b}_N = \lambda_{N-1}(1+1/1000)\f$. Finally the function trivially
        calculates the wavelength bin widths from the bin borders. */
    void setWavelengths(const Array& lambdav);

    //======================== Other Functions =======================

public:
    /** This function returns the number of bins, \f$N\f$, in the grid (or equivalently, the number
        of representative wavelengths). */
    int numWavelengths() const;

   /** This function returns the representative wavelength \f$\lambda_\ell\f$ corresponding to the
       index \f$\ell\f$. */
    double lambda(int ell) const;

    /** This function returns the width of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{b}_{\ell+1} - \lambda^\mathrm{b}_\ell\f$. */
    double dlambda(int ell) const;

    /** This function returns the left border of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{b}_\ell\f$. */
    double lambdamin(int ell) const;

    /** This function returns the right border of the wavelength bin corresponding to the index
        \f$\ell\f$, i.e. \f$\lambda^\mathrm{b}_{\ell+1}\f$. */
    double lambdamax(int ell) const;

    /** This function returns the index \f$\ell\f$ of the wavelength bin that contains the
        specified wavelength \f$\lambda\f$, i.e. for which \f$\lambda^\mathrm{b}_\ell <= \lambda <
        \lambda^\mathrm{b}_{\ell+1}\f$. As an expection, the right border of the rightmost bin is
        considered to fall inside the bin (i.e., if \f$\lambda = \lambda^\mathrm{b}_N\f$, the
        function returns \f$\ell=N-1\f$). If \f$\lambda\f$ lies outside the wavelength range of the
        grid, the value \f$\ell=-1\f$ is returned. */
    int ell(double lambda) const;

    /** This function returns the entire table of \f$N\f$ wavelength grid points. */
    const Array& lambdav() const;

    /** This function returns the entire table of \f$N\f$ wavelength bin widths. */
    const Array& dlambdav() const;

    /** This function returns the entire table of \f$N+1\f$ wavelength bin border points. */
    const Array& blambdav() const;

    //======================== Data Members ========================

private:
    // subclasses should fill these tables by calling setWavelengths() in their setupSelfBefore()
    Array _lambdav;     // N representative wavelengths
    Array _dlambdav;    // N wavelength bin widths
    Array _blambdav;    // N+1 wavelength bin border points
};

//////////////////////////////////////////////////////////////////////

#endif

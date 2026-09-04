/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GAUSSIANLINESSED_HPP
#define GAUSSIANLINESSED_HPP

#include "Array.hpp"
#include "TabulatedSED.hpp"

////////////////////////////////////////////////////////////////////

/** GaussianLinesSED is an abstract class for representing spectral energy distributions that
    consist of one or more Gaussian emission lines. Each line \f$i\f$ is centered on a wavelength
    \f$\lambda_i\f$ and has a width set by a velocity dispersion \f$s_i\f$ configured by the user
    for that line, reflecting the thermal (or other) sub-grid motion of the emitting material, in
    exactly the same way as the LyaGaussianSED class handles a single such line. Using the photon
    velocity shift \f[ v = \frac{\lambda-\lambda_i}{\lambda_i}\,c \f] as the spectral variable,
    line \f$i\f$ contributes a specific luminosity \f[ L_{\lambda,i}(\lambda) =
    \frac{L_i}{s_i\,\sqrt{2\pi}} \exp\left(-\frac{v^2}{2s_i^2}\right) \f] where \f$L_i\f$ is the
    (relative) luminosity of the line. The %SED represented by this class is the sum of the
    contributions of all configured lines, \f$L_\lambda(\lambda) = \sum_i L_{\lambda,i}(\lambda)\f$,
    properly handling any overlap between neighbouring lines.

    The subclass must provide the center wavelength, the velocity dispersion, and the relative
    luminosity for each line, in any order -- this class sorts the lines by wavelength before
    further processing, so the lines may be listed in any convenient order in the input. This
    abstract class handles everything else.

    Because a set of lines can be spread over an arbitrarily large wavelength range while each
    individual line is often extremely narrow relative to that range, this class cannot simply
    tabulate the summed profile on a single regular wavelength grid covering the full range: this
    would require an unreasonably large number of grid points to properly resolve the narrowest
    line, or would grossly undersample that line if the number of points is kept reasonable.
    Instead, for each line \f$i\f$ this class builds a private linear grid limited to the range
    \f$\lambda_i(1\pm9s_i/c)\f$ -- beyond which the Gaussian contributes negligibly; see the
    LyaGaussianSED class for a discussion of this limit -- with a resolution proportional to that
    line's own dispersion. The %SED is then tabulated on the union of all of these private grids
    (merging grid points that end up extremely close together), evaluating the sum of all Gaussians
    at each grid point. This automatically leads to a higher resolution wherever lines overlap, and
    to just the two edge points of the nearest line(s) across the (possibly wide) empty ranges in
    between non-overlapping lines or line clusters, without requiring any special-case logic. */
class GaussianLinesSED : public TabulatedSED
{
    ITEM_ABSTRACT(GaussianLinesSED, TabulatedSED,
                  "a spectral energy distribution consisting of one or more Gaussian emission lines")
        ATTRIBUTE_TYPE_DISPLAYED_IF(GaussianLinesSED, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function must be implemented in each subclass to return the center wavelength, the
        velocity dispersion, and the relative luminosity for each of the Gaussian lines in the
        %SED. The function must guarantee that the three arrays have the same, nonzero size. The
        order of the lines is not important because this abstract class sorts the lines by
        wavelength. Constant scaling of the luminosities is not important because the %SED will be
        normalized by this abstract class. */
    virtual void getLineProperties(Array& lambdav, Array& dispersionv, Array& Lv) const = 0;

    /** This function obtains the properties of the individual lines from the subclass and
        tabulates the sum of the corresponding Gaussian profiles on an appropriate wavelength grid,
        as described in the class header. */
    void getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BAND_HPP
#define BAND_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
#include "Range.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of a Band subclass represents the transmission curve of a particular observational
    filter as a function of wavelength. Examples include the standard Johnson filters and the
    transmission curves in each broadband for actual observatories such as GALEX, SDSS or Herschel.

    This abstract base class implements the key operations offered by all Band objects: obtaining
    the transmission at a given wavelength, calculating the mean specific luminosity for a given
    %SED after convolving it with the transmission curve, and determining the pivot wavelength (a
    characteristic wavelength of the band at which the mean specific luminosity can be easily
    converted between wavelength and frequency representations). Subclasses are responsible only
    for acquiring and serving the data points defining the transmission curve for the band.

    Refer to the appendix in Camps et al. 2016 (MNRAS 462, 1057-1075) for a brief introduction of
    the relevant concepts and a derivation of the corresponding formulas. For the purposes of this
    class, a band is defined through its transmission curve \f$T(\lambda)\f$. The normalized
    transmission (or the transmission probablity) at a given wavelength can then be written as \f[
    T_\text{norm}(\lambda) = \frac{ T(\lambda) } { \int T(\lambda) \,\mathrm{d}\lambda }. \f]

    Given a spectral energy distribution \f$L_\lambda(\lambda)\f$, the mean specific luminosity
    \f$\left<L_\lambda\right>\f$ can then be obtained through \f[ \left<L_\lambda\right> = \frac{
    \int L_\lambda(\lambda)T(\lambda) \,\mathrm{d}\lambda } { \int T(\lambda) \,\mathrm{d}\lambda
    }. \f]

    Furthermore, the pivot wavelength is given by \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int
    T(\lambda) \,\mathrm{d}\lambda } { \int T(\lambda) \,\mathrm{d}\lambda/\lambda^2 } }. \f]

    As set forth by Camps et al. 2016, for energy measuring devices (bolometers) the total system
    transmission \f$T(\lambda)\f$ is usually given by the instrument designers, and it can be used
    directly for this class. For photon counters (including all instruments in the UV, optical and
    near infrared), the total system response \f$R(\lambda)\f$ is usually given instead. The
    transmission curve needed for this class can be derived from the response curve through \f[
    T(\lambda) = \lambda\,R(\lambda). \f] This operation must be performed in a preprocessing step
    before handing the data to this class.

    Also, to further simplify the above formulas, this class requires that the transmission curve
    provided by a subclass is normalized to unity, i.e. \f[ \int T(\lambda) \,\mathrm{d}\lambda =
    1. \f] */
class Band : public SimulationItem
{
    ITEM_ABSTRACT(Band, SimulationItem, "a wavelength band (transmission curve)")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains pointers to the data representing the transmission curve from the
        subclass and precomputes various values for later use. */
    void setupSelfAfter() override;

    //============== Functions implemented here ============

public:
    /** This function returns the wavelength range in which the transmission for this band may be
        nonzero. */
    Range wavelengthRange() const;

    /** This function returns the normalized transmission \f$T_\text{norm}(\lambda)\f$ for this
        band at the specified wavelength. See the class header for the relevant formulas. */
    double transmission(double wavelength) const;

    /** This function returns the mean specific luminosity \f$\left<L_\lambda\right>\f$ for a given
        %SED after convolving it with the transmission curve for this band. See the class header
        for the relevant formulas. */
    double meanSpecificLuminosity(const Array& lambdav, const Array& pv) const;

    /** This function returns the pivot wavelength for this band, i.e. the wavelength at which the
        mean specific luminosity can be easily converted between wavelength and frequency
        representations. See the class header for the relevant formulas. */
    double pivotWavelength() const;

    //============== Functions to be implemented in subclasses ============

protected:
    /** This function returns the number of elements in the wavelength and transmission data arrays
        held by a subclass. */
    virtual size_t dataSize() const = 0;

    /** This function returns a pointer to the first wavelength in the corresponding array held by
        a subclass. The number of elements in this array can be obtained through the dataSize()
        function. */
    virtual const double* wavelengthData() const = 0;

    /** This function returns a pointer to the first transmission value in the corresponding array
        held by a subclass. The number of elements in this array can be obtained through the
        dataSize() function. The transmission values must be normalized as described in the class
        header. */
    virtual const double* transmissionData() const = 0;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    size_t _size{0};
    const double* _lambdav{nullptr};
    const double* _transv{nullptr};
    double _pivot{0.};
};

////////////////////////////////////////////////////////////////////

#endif

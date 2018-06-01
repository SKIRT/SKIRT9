/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SED_HPP
#define SED_HPP

#include "SimulationItem.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** An instance of a SED subclass represents a spectral energy distribution \f$L_\lambda\f$, i.e.
    power per unit of wavelength. By definition, the distribution is normalized to unity, i.e.
    integrating over all wavelengths yields one. There are two key member functions that each SED
    subclass should have: a function returning the specific luminosity at a given wavelength, and a
    function drawing a random wavelength from the spectral energy distribution.

    This abstract base class just defines an interface that must be implemented by each subclass.
    */
class SED : public SimulationItem
{
    ITEM_ABSTRACT(SED, SimulationItem, "a spectral energy distribution")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets the wavelength range for the %SED. It must be called by its owning
        source during setup \em before the setupSelfBefore() function is invoked. */
    void setWavelengthRange(double minWavelength, double maxWavelength);

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    /** This function returns the minimum of the %SED wavelength range. */
    double minWavelength() const { return _minWavelength; }

    /** This function returns the maximum of the %SED wavelength range. */
    double maxWavelength() const { return _maxWavelength; }

    //======================== Other Functions =======================

public:
    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. power per
        unit of wavelength) at the specified wavelength, or zero if the wavelength is outside of
        the distribution's spectral range. */
    virtual double specificLuminosity(double wavelength) const = 0;

    /** This function draws a random wavelength from the normalized spectral energy distribution
        represented by this object. */
    virtual double generateWavelength() const = 0;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data members initialized during setup by setWavelengthRange()
    double _minWavelength{0};
    double _maxWavelength{0};

    // data member initialized during setup
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif

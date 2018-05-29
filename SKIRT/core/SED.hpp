/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SED_HPP
#define SED_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of a SED subclass represents a spectral energy distribution \f$L_\lambda\f$, i.e.
    power per unit of wavelength. By definition, the distribution is normalized to unity, i.e.
    integrating over all wavelengths yields one.

    This abstract base class just defines an interface that must be implemented by each subclass.
    */
class SED : public SimulationItem
{
    ITEM_ABSTRACT(SED, SimulationItem, "a spectral energy distribution")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the normalized specific luminosity \f$L_\lambda\f$ (i.e. power per
        unit of wavelength) at the specified wavelength, or zero if the wavelength is outside of
        the distribution's spectral range. */
    double specificLuminosity(double wavelength) const;
};

////////////////////////////////////////////////////////////////////

#endif

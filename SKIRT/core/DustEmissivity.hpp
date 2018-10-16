/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIVITY_HPP
#define DUSTEMISSIVITY_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
class MaterialMix;

//////////////////////////////////////////////////////////////////////

/** The DustEmissivity class is the abstract base class for objects that calculate the wavelength
    dependent emissivity spectrum of a particular dust mix in a given radiation field. The
    emissivity can be determined from the optical properties of the dust mixture and the
    interstellar radiation field (both specified as arguments by the caller), and some additional
    assumptions. DustEmissivity subclasses implement various assumptions. */
class DustEmissivity : public SimulationItem
{
    ITEM_ABSTRACT(DustEmissivity, SimulationItem, "a dust emissivity calculator")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dust emissivity spectrum \f$\varepsilon_{\ell'}\f$ for the
        specified dust mix residing in a radiation field with the specified mean intensities
        \f$J_\ell\f$. The input radation field is assumed to be discretized on the radiation field
        wavelength grid as returned by the Configuration::radiationFieldWLG() function. The output
        emission spectrum is discretized on the emission spectrum wavelength grid as returned by
        the Configuration::emissionSpectrumWLG() function. */
    virtual Array emissivity(const MaterialMix* mix, const Array& Jv) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

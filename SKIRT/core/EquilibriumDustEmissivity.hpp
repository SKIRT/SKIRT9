/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EQUILIBRIUMDUSTEMISSIVITY_HPP
#define EQUILIBRIUMDUSTEMISSIVITY_HPP

#include "DustEmissivity.hpp"

//////////////////////////////////////////////////////////////////////

/** The EquilibriumDustEmissivity class calculates the emissivity of a particular dust mix in a given
    radiation field, assuming that the dust grains are in local thermal equilibrium. Under this
    assumption (which is valid for large grains) the dust emits as a modified blackbody, with a
    different equilibrium temperature for every population in the mixture if it is a multi-grain
    mixture. The emissivity in an interstellar radiation field \f$J_\lambda\f$ can then immediately
    be calculated as \f[ \varepsilon_\lambda = \frac{1}{\mu} \sum_{c=0}^{N_{\text{pop}}-1}
    \varsigma_{\lambda,c}^{\text{abs}}\, B_\lambda(T_c) \f] with \f$\mu\f$ the total dust mass of
    the dust mix, \f$\varsigma_{\lambda,c}^{\text{abs}}\f$ the absorption cross section of the
    \f$c\f$'th dust population, and \f$T_c\f$ the equilibrium temperature of that population,
    defined by the balance equation \f[ \int_0^\infty \varsigma_{\lambda,c}^{\text{abs}}\,
    J_\lambda\, {\text{d}}\lambda = \int_0^\infty \varsigma_{\lambda,c}^{\text{abs}}\,
    B_\lambda(T_c)\, {\text{d}}\lambda. \f] */
class EquilibriumDustEmissivity : public DustEmissivity
{
    ITEM_CONCRETE(EquilibriumDustEmissivity, DustEmissivity,
                  "Local thermal equilibrium (LTE) dust emissivity calculator")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dust emissivity spectrum \f$\varepsilon_{\ell'}\f$ for the
        specified dust mix residing in a radiation field with the specified mean intensities
        \f$J_\ell\f$. The input and output arrays are discretized on the wavelength grids returned
        by the Configuration::radiationFieldWLG() and Configuration::emissionSpectrumWLG()
        functions, repectively. */
    Array emissivity(const MaterialMix* mix, const Array& Jv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EQUILIBRIUMDUSTEMISSIVITY_HPP
#define EQUILIBRIUMDUSTEMISSIVITY_HPP

#include "DustEmissivity.hpp"

//////////////////////////////////////////////////////////////////////

/** The EquilibriumDustEmissivity class calculates the emissivity of a particular dust mix in a
    given radiation field, assuming that the dust grains are in local thermal equilibrium. Under
    this assumption (which is usually only valid for large grains) the dust emits as a modified
    blackbody, with an equilibrium temperature depending on material composition and grain size.

    For dust mixes described by a single representative dust grain, the class uses the properties
    of that grain to obtain a single equilibrium temperature and resulting emissivity. For dust
    mixes described by multiple representative grains (for various material compositions and grain
    size bins), the class accumulates the emissivities at the equilibrium temperature for each of
    representative grains, resulting in a more accurate result (still limited by the assumption of
    local thermal equilibrium). Expert users can force the class to use a single representative
    grain even for dust mixes described by multiple grains.

    In formula form, the emissivity of a dust mix in an embedding radiation field \f$J_\lambda\f$
    can be calculated as \f[ \varepsilon_\lambda = \frac{1}{\mu} \sum_{b=0}^{N_{\text{bins}}-1}
    \varsigma_{\lambda,b}^{\text{abs}}\, B_\lambda(T_b) \f] with \f$\mu\f$ the total dust mass of
    the dust mix, \f$\varsigma_{\lambda,b}^{\text{abs}}\f$ the absorption cross section of the
    \f$b\f$'th representative grain, and \f$T_b\f$ the equilibrium temperature of that grain,
    defined by the balance equation \f[ \int_0^\infty \varsigma_{\lambda,b}^{\text{abs}}\,
    J_\lambda\, {\text{d}}\lambda = \int_0^\infty \varsigma_{\lambda,b}^{\text{abs}}\,
    B_\lambda(T_b)\, {\text{d}}\lambda. \f] */
class EquilibriumDustEmissivity : public DustEmissivity
{
    ITEM_CONCRETE(EquilibriumDustEmissivity, DustEmissivity,
                  "Local thermal equilibrium (LTE) dust emissivity calculator")

    PROPERTY_BOOL(useSingleGrain, "use a single representative grain even for multi-grain dust mixes")
        ATTRIBUTE_DEFAULT_VALUE(useSingleGrain, "false")
        ATTRIBUTE_DISPLAYED_IF(useSingleGrain, "Level3")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the dust emissivity spectrum \f$\varepsilon_{\ell'}\f$ for the
        specified dust mix residing in a radiation field with the specified mean intensities
        \f$J_\ell\f$. The input and output arrays are discretized on the wavelength grids returned
        by the Configuration::radiationFieldWLG() and Configuration::dustEmissionWLG() functions,
        repectively. */
    Array emissivity(const MaterialMix* mix, const Array& Jv) const override;
};

////////////////////////////////////////////////////////////////////

#endif

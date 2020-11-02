/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEANTRUSTBENCHMARKDUSTMIX_HPP
#define MEANTRUSTBENCHMARKDUSTMIX_HPP

#include "SingleGrainDustMix.hpp"

////////////////////////////////////////////////////////////////////

/** The MeanTrustBenchmarkDustMix class represents a population of identical dust grains with
    properties approximating those of a mixture of bare (i.e. non-composite) graphite, silicate and
    PAH dust grains. The size distribution of each of these dust grain populations is fine-tuned in
    such a way that the global dust properties accurately reproduce the extinction, emission and
    abundance constraints on the Milky Way. The size distributions are taken from Zubko, Dwek &
    Arendt (2004, ApJS, 152, 211) and correspond to model BARE_GR_S.

    This dust mix has been prepared by Karl Misselt for the TRUST benchmark simulations; see Camps
    et al. 2015 (A&A 580:A87) and Gordon et al. 2017 (A&A 603:A114). The data can be downloaded
    from http://www.shg.ugent.be/html/_downloads.html or http://ipag.osug.fr/RT13/RTTRUST/opa.php
    */
class MeanTrustBenchmarkDustMix : public SingleGrainDustMix
{
    /** The enumeration type indicating the scattering mode. */
    ENUM_DEF(ScatteringType, HenyeyGreenstein, MaterialPhaseFunction, SphericalPolarization)
        ENUM_VAL(ScatteringType, HenyeyGreenstein, "use the Henyey-Greenstein phase function (unpolarized)")
        ENUM_VAL(ScatteringType, MaterialPhaseFunction,
                 "use the phase function derived from actual material properties (unpolarized)")
        ENUM_VAL(ScatteringType, SphericalPolarization, "support polarization through scattering by spherical grains")
    ENUM_END()

    ITEM_CONCRETE(MeanTrustBenchmarkDustMix, SingleGrainDustMix,
                  "a TRUST benchmark dust mix (mean properties, optionally with polarization)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MeanTrustBenchmarkDustMix, "Level2")

        PROPERTY_ENUM(scatteringType, ScatteringType, "the type of scattering to be implemented")
        ATTRIBUTE_DEFAULT_VALUE(scatteringType, "HenyeyGreenstein")

    ITEM_END()

    //======================== Capabilities =======================

public:
    /** This function returns the scattering mode supported by this material mix as configured by
        the user through the scatteringType property. */
    ScatteringMode scatteringMode() const override;

    /** This function returns a flag indicating whether the material mix supports polarization
        during scattering events or not. For this dust mix, the function returns true if the user
        configured the SphericalPolarization scattering type, and false otherwise. */
    bool hasPolarizedScattering() const override;

    //======================== Other Functions =======================

protected:
    /** This function returns the name of the stored table resource tabulating the basic optical
        properties for this dust mix. */
    string resourceNameForOpticalProps() const override;

    /** This function returns the name of the stored table resource tabulating the Mueller matrix
        elements for this dust mix. */
    string resourceNameForMuellerMatrix() const override;
};

////////////////////////////////////////////////////////////////////

#endif

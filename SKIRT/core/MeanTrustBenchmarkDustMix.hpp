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
    ITEM_CONCRETE(MeanTrustBenchmarkDustMix, SingleGrainDustMix,
                  "a dust mix from the TRUST benchmark (mean properties, optionally with polarization)")

    PROPERTY_BOOL(includePolarization, "include support for polarization by scattering")
        ATTRIBUTE_DEFAULT_VALUE(includePolarization, "false")

    ITEM_END()

    //======================== Other Functions =======================

    /** This function returns the name of the stored table resource tabulating the basic optical
        properties for this dust mix. */
    string resourceNameForOpticalProps() const override;

    /** If the \em includePolarization flag is enabled, this function returns the name of the
        stored table resource tabulating the Mueller matrix elements for this dust mix. Otherwise
        it returns the empty string. */
    string resourceNameForMuellerMatrix() const override;
};

////////////////////////////////////////////////////////////////////

#endif

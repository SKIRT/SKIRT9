/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LUMINOSITYNORMALIZATION_HPP
#define LUMINOSITYNORMALIZATION_HPP

#include "SimulationItem.hpp"
class SED;

////////////////////////////////////////////////////////////////////

/** An instance of a LuminosityNormalization subclass describes the normalization of the luminosity
    for a primary source. In other words, it provides the mechanism for determining the appropriate
    scaling factor given a spectral distribution that is normalized to unity. This abstract base
    class just defines an interface that must be implemented by each subclass. */
class LuminosityNormalization : public SimulationItem
{
    ITEM_ABSTRACT(LuminosityNormalization, SimulationItem, "the luminosity normalization for a primary source")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the bolometric luminosity of a source with the normalized spectral
        distribution described by the specified SED object. */
    virtual double luminosity(SED* sed) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

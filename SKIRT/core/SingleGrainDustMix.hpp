/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEGRAINDUSTMIX_HPP
#define SINGLEGRAINDUSTMIX_HPP

#include "MaterialMix.hpp"
#include "DustExtinctionMixInterface.hpp"
#include "ScatteringMixInterface.hpp"

////////////////////////////////////////////////////////////////////

/** SingleGrainDustMix is an abstract class describing a dust mix described by a single
    representative grain, with or without support for polarization by scattering. This base class
    includes the implementations of the required functions. Subclasses must merely provide the
    names of the relevant resource files. */
class SingleGrainDustMix : public MaterialMix, public DustExtinctionMixInterface, public ScatteringMixInterface
{
    ITEM_ABSTRACT(SingleGrainDustMix, MaterialMix, "a dust mix described by a single representative grain")
    ITEM_END()

    //============= Implementing ElectronExtinctionMixInterface =============


    //============= Implementing ScatteringMixInterface =============

    /** This function returns true if this dust mix supports polarization; false otherwise. */
    bool hasScatteringPolarization() const override;
};

////////////////////////////////////////////////////////////////////

#endif

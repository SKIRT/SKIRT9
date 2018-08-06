/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONMIX_HPP
#define ELECTRONMIX_HPP

#include "MaterialMix.hpp"
#include "ElectronExtinctionMixInterface.hpp"
#include "ScatteringMixInterface.hpp"

////////////////////////////////////////////////////////////////////

/** The ElectronMix class describes the material properties for a population of electrons. */
class ElectronMix : public MaterialMix, public ElectronExtinctionMixInterface, public ScatteringMixInterface
{
    // TODO: make concrete
    ITEM_ABSTRACT(ElectronMix, MaterialMix, "a population of electrons")
    ITEM_END()

    //============= Implementing ElectronExtinctionMixInterface =============


    //============= Implementing ScatteringMixInterface =============

    /** This function returns true because the electron mix supports polarization. */
    bool hasScatteringPolarization() const override;
};

////////////////////////////////////////////////////////////////////

#endif

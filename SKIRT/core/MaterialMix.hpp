/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALMIX_HPP
#define MATERIALMIX_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** MaterialMix is the abstract base class for all classes representing the concrete material
    properties of a specific transfer medium. The MaterialMix class hierarchy supports three types
    of material: dust, electrons, and hydrogen-dominated gas. Many MaterialMix subclasses represent
    just one of these material types, however some may represent a combination of two or even all
    three of these material types.

    The MaterialMix base class offers no properties nor functions (other than those common to all
    simulation items). Subclasses are expected to implement the relevant "XxxMixInterface"
    interfaces depending on the represented material type(s) and depending on the level of support
    offered for each. For example, a dust mix describing a single representative grain may
    implement just the DustExtinctionMixInterface and the ScatteringMixInterface, while other dust
    mixes may also implement more advanced interfaces that provide properties needed to calculate
    stochastic dust heating.

    In any case, each subclass must implement at least one XxxExtinctionMixInterface where Xxx is
    Dust, Electron, or Gas, because that is how the Medium object to which the MaterialMix is
    attached determines which material type(s) is/are being described.

    It is worth noting at this point that the physical quantities used to represent the amount of
    material in a given location or in a spatial region depend on the material type as indicated in
    the table below. These differences are reflected in the corresponding interfaces.

    Material type | amount per volume at given location | amount in given region
    --------------|-------------------------------------|-----------------------
    Dust          | dust mass density (kg/m3)           | dust mass (kg)
    Electron      | electron number density (\#/m3)     | electron number (\#)
    Gas           | proton number density (\#/m3)       | proton number (\#)

    */
class MaterialMix : public SimulationItem
{
    ITEM_ABSTRACT(MaterialMix, SimulationItem, "a material mix")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

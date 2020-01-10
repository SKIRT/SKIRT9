/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALNORMALIZATION_HPP
#define MATERIALNORMALIZATION_HPP

#include "SimulationItem.hpp"
class Geometry;
class MaterialMix;

//////////////////////////////////////////////////////////////////////

/** MaterialNormalization is an abstract base class for classes that allow specifying the amount
    of material in a GeometricMedium instance in various ways. */
class MaterialNormalization : public SimulationItem
{
    ITEM_ABSTRACT(MaterialNormalization, SimulationItem, "a normalization for the amount of material")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the total number of entities and total mass in the medium, in that
        order, given a geometry and material mix in addition to user configuration options offered
        by the subclass implementing the function. */
    virtual std::pair<double, double> numberAndMass(const Geometry* geom, const MaterialMix* mix) const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif

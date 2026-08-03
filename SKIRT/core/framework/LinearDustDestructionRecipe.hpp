/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LINEARDUSTDESTRUCTIONRECIPE_HPP
#define LINEARDUSTDESTRUCTIONRECIPE_HPP

#include "DustDestructionRecipe.hpp"

////////////////////////////////////////////////////////////////////

/** LinearDustDestructionRecipe derives from DustDestructionRecipe to implement a basic dust
    destruction recipe that depends solely on the equilibrium temperature \f$T_\mathrm{eq}\f$ of
    the grain and on two cutoff temperatures \f$T_\mathrm{min}\f$ and \f$T_\mathrm{max}\f$. These
    cutoff temperatures can be configured separately for each type of grain material (silicate or
    graphite).

    Specifically, there is no destruction if \f$T_\mathrm{eq}<=T_\mathrm{min}\f$; all grains are
    destroyed if \f$T_\mathrm{eq}>=T_\mathrm{min}\f$; and in between the destruction fraction is
    obtained by linear interpolation. */
class LinearDustDestructionRecipe : public DustDestructionRecipe
{
    ITEM_CONCRETE(LinearDustDestructionRecipe, DustDestructionRecipe,
                  "a dust destruction recipe using a linear temperature dependence")

        PROPERTY_DOUBLE(minSilicateTemperature, "the temperature below which silicate grains are not destroyed")
        ATTRIBUTE_QUANTITY(minSilicateTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(minSilicateTemperature, "[3 K")
        ATTRIBUTE_MAX_VALUE(minSilicateTemperature, "30000 K]")
        ATTRIBUTE_DEFAULT_VALUE(minSilicateTemperature, "1200 K")

        PROPERTY_DOUBLE(maxSilicateTemperature, "the temperature above which all silicate grains are destroyed")
        ATTRIBUTE_QUANTITY(maxSilicateTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(maxSilicateTemperature, "[3 K")
        ATTRIBUTE_MAX_VALUE(maxSilicateTemperature, "30000 K]")
        ATTRIBUTE_DEFAULT_VALUE(maxSilicateTemperature, "1200 K")

        PROPERTY_DOUBLE(minGraphiteTemperature, "the temperature below which graphite grains are not destroyed")
        ATTRIBUTE_QUANTITY(minGraphiteTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(minGraphiteTemperature, "[3 K")
        ATTRIBUTE_MAX_VALUE(minGraphiteTemperature, "30000 K]")
        ATTRIBUTE_DEFAULT_VALUE(minGraphiteTemperature, "2000 K")

        PROPERTY_DOUBLE(maxGraphiteTemperature, "the temperature above which all graphite grains are destroyed")
        ATTRIBUTE_QUANTITY(maxGraphiteTemperature, "temperature")
        ATTRIBUTE_MIN_VALUE(maxGraphiteTemperature, "[3 K")
        ATTRIBUTE_MAX_VALUE(maxGraphiteTemperature, "30000 K]")
        ATTRIBUTE_DEFAULT_VALUE(maxGraphiteTemperature, "2000 K")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the maximum temperatures are not below the minimum
        temperatures. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the non-destroyed density fraction for a grain population with the
        specified type (graphite or silicate) and equilibrium temperature using the linear
        dependence described in the class header. The grain radius and the radiation field are not
        used. */
    double densityFraction(bool graphite, double a, const Array& Jv, double T) const override;
};

////////////////////////////////////////////////////////////////////

#endif

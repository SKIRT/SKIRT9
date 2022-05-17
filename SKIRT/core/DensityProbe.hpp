/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DENSITYPROBE_HPP
#define DENSITYPROBE_HPP

#include "SpatialGridWhenFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** DensityProbe probes the density of the medium as discretized on the spatial grid of the
    simulation. For dust components, the probe outputs mass density, while for electron and gas
    components it outputs number density. The probe can be used with any Form subclass. When
    associated with a form that projects the quantity along a path, the probe outputs mass surface
    density or number surface density (column density) depending on the material type.

    The user can select the aggregation level, i.e. whether to produce an output file per medium
    component or per medium type (dust, electrons, gas). There is also an option to decide whether
    the probe should be performed after setup or after the full simulation run. The latter option
    is meaningful if the density of the media may change during the simulation. */
class DensityProbe : public SpatialGridWhenFormProbe
{
    /** The enumeration type indicating how to aggregate the output: per medium component or per
        medium type (dust, electrons, gas). */
    ENUM_DEF(Aggregation, Component, Type)
        ENUM_VAL(Aggregation, Component, "per medium component")
        ENUM_VAL(Aggregation, Type, "per medium type (dust, electrons, gas)")
    ENUM_END()

    ITEM_CONCRETE(DensityProbe, SpatialGridWhenFormProbe, "internal spatial grid: density of the medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DensityProbe, "Medium&SpatialGrid")

        PROPERTY_ENUM(aggregation, Aggregation, "how to aggregate the density")
        ATTRIBUTE_DEFAULT_VALUE(aggregation, "Type")
        ATTRIBUTE_DISPLAYED_IF(aggregation, "Level2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

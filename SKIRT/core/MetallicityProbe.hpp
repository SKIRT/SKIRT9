/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef METALLICITYPROBE_HPP
#define METALLICITYPROBE_HPP

#include "SpatialGridWhenFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** MetallicityProbe probes the metallicity of gas medium components as discretized on the spatial
    grid of the simulation. It works only for gas components that store a metallicity value in the
    medium state. Electrons do not have metallicity. Dust components do not store metallicity in
    the medium state even if the value is imported from an input file; the metallicity value is
    discarded immediately after being used for calculating the dust density.

    A seperate output file is produced for each gas medium component that stores a metallicity
    value in the medium state. The probe can be used with any Form subclass. When associated with a
    form that projects the quantity along a path, the value is density-weighted. */
class MetallicityProbe : public SpatialGridWhenFormProbe
{
    ITEM_CONCRETE(MetallicityProbe, SpatialGridWhenFormProbe, "internal spatial grid: metallicity of the gas medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MetallicityProbe, "Medium&SpatialGrid&GasMix")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

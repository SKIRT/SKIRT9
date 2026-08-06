/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VELOCITYPROBE_HPP
#define VELOCITYPROBE_HPP

#include "SpatialGridFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** VelocityProbe probes the bulk velocity of the medium as discretized on the spatial grid of the
    simulation. The probe can be used with any Form subclass. When associated with a form that
    projects the quantity along a path, the value is density-weighted. */
class VelocityProbe : public SpatialGridFormProbe
{
    ITEM_CONCRETE(VelocityProbe, SpatialGridFormProbe, "internal spatial grid: bulk velocity of the medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(VelocityProbe, "Level2&SpatialGrid&MediumVelocity")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

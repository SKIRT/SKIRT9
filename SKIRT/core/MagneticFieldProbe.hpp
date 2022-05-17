/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MAGNETICFIELDPROBE_HPP
#define MAGNETICFIELDPROBE_HPP

#include "SpatialGridFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** MagneticFieldProbe probes the magnetic field associated with the medium as discretized on the
    spatial grid of the simulation. The probe can be used with any Form subclass. When associated
    with a form that projects the quantity along a path, the value is density-weighted. */
class MagneticFieldProbe : public SpatialGridFormProbe
{
    ITEM_CONCRETE(MagneticFieldProbe, SpatialGridFormProbe,
                  "internal spatial grid: magnetic field associated with the medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MagneticFieldProbe, "Level3&SpatialGrid&MagneticField")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

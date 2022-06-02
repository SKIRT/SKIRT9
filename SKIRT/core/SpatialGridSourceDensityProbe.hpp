/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDSOURCEDENSITYPROBE_HPP
#define SPATIALGRIDSOURCEDENSITYPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SpatialGridSourceDensityProbe is intended for specialty situations where one wants to get
    information on the discretization of the source distribution over the spatial grid in the
    simulation. This makes sense and is supported only for source components of type
    GeometricSource. In other cases, the spatial distribution of the source is either uninteresting
    (such as for background sources) or wavelength dependent (such as for imported sources). Also,
    each geometric source component is treated seperately to avoid the need for taking into account
    the luminosity normalization of the respective source components.

    When the stated requirements are met, the probe outputs a text column file, named
    <tt>prefix_probe_sourcedens.dat</tt>, which contains a line for each cell in the spatial grid.
    Each line contains the cell index and the normalized density for each GeometricSource component
    sampled at the cell center. */
class SpatialGridSourceDensityProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(SpatialGridSourceDensityProbe, SpecialtyProbe,
                  "specialty: primary source density discretized on spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpatialGridSourceDensityProbe, "Level3&GeometricSource&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

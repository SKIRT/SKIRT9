/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTDUSTTEMPERATURECUTSPROBE_HPP
#define DEFAULTDUSTTEMPERATURECUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultDustTemperatureCutsProbe outputs FITS files (named <tt>prefix_probe_dust_T_XX.fits</tt>)
    with cuts through an indicative dust temperature along the coordinate planes. The number of
    data files written depends on the geometry and material contents of the media system. For
    spherical symmetry only the intersection with the xy plane is written, for axial symmetry the
    intersections with the xy and xz planes are written, and for general geometries all three
    intersections are written. Each of the output files contains a map with 1024 x 1024 pixels, and
    covers a field of view equal to the total extension of the spatial grid in the simulation.

    The indicative dust temperature assigned to each output pixel is obtained as follows. First the
    probe determines the cell in the simulation's spatial grid "containing" the pixel according to
    its position in the cut being considered. For each material mix of type dust present in the
    cell, or if applicable, for each dust population in these mixes, the probe then calculates the
    equilibrium temperature that would be reached when the dust is embedded in the radiation field
    tracked by the simulation for the cell. This is achieved by solving the energy balance equation
    under LTE (local thermal equilibrium) assumptions. The resulting temperatures are finally
    averaged over the dust populations in each mix (weighed by the relative mass in the mix) over
    and all dust components present in the spatial cell (weighed by relative mass in the cell).

    Note that the indicative dust temperature does not correspond to a physical temperature. The
    LTE assumption is almost certainly unjustified for a relevant portion of the dust grains
    (depending on the embedding radiation field), and even when ignoring this problem, averaging
    temperatures over dust populations and dust mixes has no clear-cut physical interpretation. */
class DefaultDustTemperatureCutsProbe : public Probe
{
    ITEM_CONCRETE(DefaultDustTemperatureCutsProbe, Probe,
                  "cuts of the indicative dust temperature along the coordinate axes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultDustTemperatureCutsProbe, "Dust&SpatialGrid&RadiationField")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif

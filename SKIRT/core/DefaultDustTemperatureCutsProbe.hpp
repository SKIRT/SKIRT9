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
    its position in the cut being considered. The probe then calculates the indicative dust
    temperature for that cell by averaging the LTE equilibrium temperatures for the various dust
    mixes present in the cell. Note that the indicative dust temperature does not really correspond
    to a physical temperature. For more information about the indicative dust temperature, refer to
    the MediumSystem::indicativeDustTemperature() function. */
class DefaultDustTemperatureCutsProbe : public Probe
{
    ITEM_CONCRETE(DefaultDustTemperatureCutsProbe, Probe,
                  "cuts of the indicative dust temperature along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultDustTemperatureCutsProbe,
                                    "Dust&SpatialGrid&RadiationField&Panchromatic")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif

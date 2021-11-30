/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTELECTRONTEMPERATURECUTSPROBE_HPP
#define DEFAULTELECTRONTEMPERATURECUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultElectronTemperatureCutsProbe outputs FITS files (named
    <tt>prefix_probe_dust_T_XX.fits</tt>) with cuts through the electron temperature defined by the
    input model along the coordinate planes. The number of data files written depends on the
    geometry and material contents of the media system. For spherical symmetry only the
    intersection with the xy plane is written, for axial symmetry the intersections with the xy and
    xz planes are written, and for general geometries all three intersections are written. Each of
    the output files contains a map with 1024 x 1024 pixels, and covers a field of view equal to
    the total extent of the spatial grid in the simulation. */
class DefaultElectronTemperatureCutsProbe : public Probe
{
    ITEM_CONCRETE(DefaultElectronTemperatureCutsProbe, Probe,
                  "cuts of the electron temperature along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultElectronTemperatureCutsProbe, "Level2&ElectronMix")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing at the end of the setup phase. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTGASTEMPERATURECUTSPROBE_HPP
#define DEFAULTGASTEMPERATURECUTSPROBE_HPP

#include "StateProbe.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultGasTemperatureCutsProbe outputs FITS files (named <tt>prefix_probe_dust_T_XX.fits</tt>)
    with cuts through the indicative gas temperature along the coordinate planes. The number of
    data files written depends on the geometry and material contents of the media system. For
    spherical symmetry only the intersection with the xy plane is written, for axial symmetry the
    intersections with the xy and xz planes are written, and for general geometries all three
    intersections are written. Each of the output files contains a map with 1024 x 1024 pixels, and
    covers a field of view equal to the total extent of the spatial grid in the simulation. */
class DefaultGasTemperatureCutsProbe : public StateProbe
{
    ITEM_CONCRETE(DefaultGasTemperatureCutsProbe, StateProbe, "cuts of the gas temperature along the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultGasTemperatureCutsProbe, "Level2&Gas&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs the probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

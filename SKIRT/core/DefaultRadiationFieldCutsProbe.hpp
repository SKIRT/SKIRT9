/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTRADIATIONFIELDCUTSPROBE_HPP
#define DEFAULTRADIATIONFIELDCUTSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** DefaultRadiationFieldCutsProbe outputs FITS files (named <tt>prefix_probe_J_XX.fits</tt>) with
    cuts through the mean radiation field intensity along the coordinate planes. The number of data
    files written depends on the geometry and material contents of the media system. For spherical
    symmetry only the intersection with the xy plane is written, for axial symmetry the
    intersections with the xy and xz planes are written, and for general geometries all three
    intersections are written.

    Each of the output files actually contains a datacube with a map for each bin in the wavelength
    grid returned by the Configuration::radiationFieldWavelengthGrid() function. The maps have 1024
    x 1024 pixels, and cover a field of view equal to the total extension of the spatial grid in
    the simulation. */
class DefaultRadiationFieldCutsProbe : public Probe
{
    ITEM_CONCRETE(DefaultRadiationFieldCutsProbe, Probe,
                  "cuts of the mean radiation field intensity along the coordinate axes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DefaultRadiationFieldCutsProbe, "Medium&SpatialGrid&RadiationField")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif

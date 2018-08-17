/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDCONVERGENCEPROBE_HPP
#define SPATIALGRIDCONVERGENCEPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** SpatialGridConvergenceProbe outputs a text file named <tt>prefix_probe_convergence.dat</tt>
    with convergence information on the spatial grid for each material type in the medium system.
    The file is formatted for human consumption (not in column text format) and is intended as a
    basic sanity check on the configuration of the simulation.

    For each material type, the file lists the total mass and the optical depth along the major
    coordinate axes at a given wavelength. These numbers can be compared to the values expected for
    the model.

    Furthermore, in each case, the file provides the value as represented by the input model
    defined by the media system, and also the grid-discretized value as obtained from the
    finite-resolution spatial grid in the simulation. A comparison of both sets of values offers a
    first indication of whether the configured spatial grid properly captures the material mass in
    the simulation (in the ideal case, there would be no difference between both sets of values).

    Finally, the file includes some basic statistics on the diagonal optical depths of the spatial
    grid cells: the largest and average diagonal optical depth, and the 90% percentile. */
class SpatialGridConvergenceProbe : public Probe
{
    ITEM_CONCRETE(SpatialGridConvergenceProbe, Probe, "convergence information on the spatial grid")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to list the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 A")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif

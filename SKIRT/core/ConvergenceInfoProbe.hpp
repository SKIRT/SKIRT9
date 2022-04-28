/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONVERGENCEINFOPROBE_HPP
#define CONVERGENCEINFOPROBE_HPP

#include "SpecialtyWavelengthProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ConvergenceInfoProbe outputs a text file named <tt>prefix_probe_convergence.dat</tt>
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
class ConvergenceInfoProbe : public SpecialtyWavelengthProbe
{
    ITEM_CONCRETE(ConvergenceInfoProbe, SpecialtyWavelengthProbe, "convergence: information on the spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ConvergenceInfoProbe, "Medium&SpatialGrid")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

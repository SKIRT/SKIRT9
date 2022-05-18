/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DEFAULTCUTSFORM_H
#define DEFAULTCUTSFORM_H

#include "SpatialGridForm.hpp"

//////////////////////////////////////////////////////////////////////

/** DefaultCutsForm represents a spatial grid-specific probe form. Refer to the ProbeFormBridge
    class for more information about probes and forms.

    This particular form outputs FITS files representing cuts through the spatial domain of the
    simulated model along the coordinate planes. The field of view of each cut covers the extent of
    the spatial grid in the simulation in the relevant directions. Each cut has 1024 x 1024 pixels.
    The number of data files written depends on the geometry of the medium system. For spherical
    symmetry only the intersection with the xy plane is written, for axial symmetry the
    intersections with the xy and xz planes are written, and for general geometries all three
    intersections are written.

    Each of the FITS files contains a number of image frames corresponding to the number of
    components in a value of the quantity being probed (i.e. 1 for scalar quantities, 3 for vector
    quantities, and N for compound quantities). In case of a vector quantity, the three image
    frames representing the velocity vector components in the frame of the cut, i.e. the two
    components projected on the x and y axes of the cut and the component perpendicular to the cut,
    where positive values indicate vectors pointing away from the viewer. */
class DefaultCutsForm : public SpatialGridForm
{
    ITEM_CONCRETE(DefaultCutsForm, SpatialGridForm, "default planar cuts along the coordinate planes")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded form of
        this type (as opposed to selected through the ski file). Before the constructor returns,
        the newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit DefaultCutsForm(SimulationItem* parent);

    //============= Other functions =============

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;
};

//////////////////////////////////////////////////////////////////////

#endif

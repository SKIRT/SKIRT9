/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDPLOTPROBE_HPP
#define SPATIALGRIDPLOTPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SpatialGridPlotProbe outputs text data files that allow plotting the structure of the spatial
    grid configured for the simulation.

    The number of data files written depends on the dimension of the spatial grid: for spherical
    symmetry only the intersection with the xy plane is written, for axial symmetry the
    intersections with the xy and xz planes are written, and for general geometries all three
    intersections are written. In the latter case, an extra file with three-dimensional information
    is written as well.

    The output files are called <tt>prefix_probe_grid_XXX.dat</tt>, where XXX is replaced by "xy",
    "xz", "yz" or "xyz" depending on the file under consideration. Within a file, each line
    contains two coordinates seperated by whitespace or is empty. Consecutive nonempty lines
    represent a sequence of "lineto" commands; an empty line marks a "moveto" command. */
class SpatialGridPlotProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(SpatialGridPlotProbe, SpecialtyProbe, "properties: data files for plotting the structure of the grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpatialGridPlotProbe, "SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

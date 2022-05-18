/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALCELLPROPERTIESPROBE_HPP
#define SPATIALCELLPROPERTIESPROBE_HPP

#include "SpecialtyWavelengthProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SpatialCellPropertiesProbe outputs a column text file, named
    <tt>prefix_probe_cellprops.dat</tt>, which contains a line for each cell in the spatial grid of
    the simulation. Each line contains columns representing the following cell properties:
    cell index, x,y,z coordinates of the cell center, volume, total optical depth of the cell
    diagonal at a user-configured wavelength (for all material types combined), dust mass density,
    electron number density, and (gas) hydrogen number density. */
class SpatialCellPropertiesProbe : public SpecialtyWavelengthProbe
{
    ITEM_CONCRETE(SpatialCellPropertiesProbe, SpecialtyWavelengthProbe, "properties: basic info for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpatialCellPropertiesProbe, "Medium&SpatialGrid")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALCELLPROPERTIESPROBE_HPP
#define SPATIALCELLPROPERTIESPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** SpatialCellPropertiesProbe outputs a column text file, named
    <tt>prefix_probe_cellprops.dat</tt>, which contains a line for each cell in the spatial grid of
    the simulation. Each line contains columns representing the following cell properties:
    cell index, x,y,z coordinates of the cell center, volume, total optical depth of the cell
    diagonal at a user-configured wavelength (for all material types combined), dust mass density,
    electron number density, and (gas) hydrogen number density. */
class SpatialCellPropertiesProbe : public Probe
{
    ITEM_CONCRETE(SpatialCellPropertiesProbe, Probe, "relevant properties for all spatial cells")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpatialCellPropertiesProbe, "Level2&SpatialGrid")

    PROPERTY_DOUBLE(wavelength, "the wavelength at which to list the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 Angstrom")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif

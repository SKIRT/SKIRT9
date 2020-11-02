/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALCELLPROPERTIESPROBE_HPP
#define SPATIALCELLPROPERTIESPROBE_HPP

#include "AbstractWavelengthProbe.hpp"

////////////////////////////////////////////////////////////////////

/** SpatialCellPropertiesProbe outputs a column text file, named
    <tt>prefix_probe_cellprops.dat</tt>, which contains a line for each cell in the spatial grid of
    the simulation. Each line contains columns representing the following cell properties:
    cell index, x,y,z coordinates of the cell center, volume, total optical depth of the cell
    diagonal at a user-configured wavelength (for all material types combined), dust mass density,
    electron number density, and (gas) hydrogen number density. */
class SpatialCellPropertiesProbe : public AbstractWavelengthProbe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(SpatialCellPropertiesProbe, AbstractWavelengthProbe, "relevant properties for all spatial cells")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpatialCellPropertiesProbe, "Level2&Medium&SpatialGrid")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "when to probe the medium state")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "HasDynamicState")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. It produces output only if the \em
        probeAfter property is set to Setup. */
    void probeSetup() override;

    /** This function performs probing after all photon packets have been emitted and detected. It
        produces output only if the \em probeAfter property is set to Run. */
    void probeRun() override;

private:
    /** This function performs the probing; it is called from probeSetup() or probeRun() depending
        on the value of the \em probeAfter property. */
    void probe();
};

////////////////////////////////////////////////////////////////////

#endif

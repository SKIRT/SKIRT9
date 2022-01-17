/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIATIONFIELDPERCELLPROBE_HPP
#define RADIATIONFIELDPERCELLPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** RadiationFieldPerCellProbe outputs a column text file (named <tt>prefix_probe_J.dat</tt>)
    listing the mean radiation field intensity for each cell in the spatial grid of the simulation.

    Specifically, the output file contains a line for each cell in the spatial grid of the
    simulation. The first columns specifies the cell index, and subsequent columns list the mean
    radiation field intensity for each bin in the wavelength grid returned by the
    Configuration::radiationFieldWLG() function.

    The probe offers an option to output a separate text column file with details on the radiation
    field wavelength grid. For each wavelength bin, the file lists the characteristic wavelength,
    the wavelength bin width, and the left and right borders of the bin. */
class RadiationFieldPerCellProbe : public Probe
{
    ITEM_CONCRETE(RadiationFieldPerCellProbe, Probe, "the mean radiation field intensity for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(RadiationFieldPerCellProbe, "Level2&SpatialGrid&RadiationField")

        PROPERTY_BOOL(writeWavelengthGrid, "output a text file with the radiation field wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(writeWavelengthGrid, "false")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif

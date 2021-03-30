/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIATIONFIELDATPOSITIONSPROBE_HPP
#define RADIATIONFIELDATPOSITIONSPROBE_HPP

#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** RadiationFieldAtPositionsProbe outputs a column text file (named <tt>prefix_probe_J.dat</tt>)
    listing the mean radiation field intensity at the positions listed in a user-provided input
    file. This allows, for example, to probe the radiation field at the positions listed in a
    Particle or VoronoiMesh import file used for defining sources or media by configuring this
    probe with the same import file (any columns other than the position coordinates will be
    ignored).

    The probe expects three columns in the input column text file, respectively representing the x,
    y and z coordinates of a position in the spatial domain of the simulation. In case the input
    file has no unit specifications, the default units are parsec. Refer to the description of the
    TextInFile class for more information on overall formatting and on how to include header lines
    specifying the units for each column.

    The output file contains a line for each (non-empty and non-comment) line in the input file, in
    the same order. The first three columns repeat the x, y and z coordinates of the position (in
    the configured output units, which may differ from the input units). Subsequent columns list
    the mean radiation field intensity for each bin in the wavelength grid returned by the
    Configuration::radiationFieldWLG() function. The radiation field at positions outside of
    simulation's spatial grid is assumed to be zero.

    The probe offers an option to output a separate text column file with details on the radiation
    field wavelength grid. For each wavelength bin, the file lists the characteristic wavelength,
    the wavelength bin width, and the left and right borders of the bin. */
class RadiationFieldAtPositionsProbe : public Probe
{
    ITEM_CONCRETE(RadiationFieldAtPositionsProbe, Probe, "the mean radiation field intensity at imported positions")
        ATTRIBUTE_TYPE_DISPLAYED_IF(RadiationFieldAtPositionsProbe, "Level2&Medium&SpatialGrid&RadiationField")

        PROPERTY_STRING(filename, "the name of the file listing the positions")

        PROPERTY_STRING(useColumns, "a list of names corresponding to columns in the file to be imported")
        ATTRIBUTE_DEFAULT_VALUE(useColumns, "")
        ATTRIBUTE_REQUIRED_IF(useColumns, "false")
        ATTRIBUTE_DISPLAYED_IF(useColumns, "Level3")

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

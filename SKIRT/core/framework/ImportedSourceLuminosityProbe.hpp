/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDSOURCELUMINOSITYPROBE_HPP
#define IMPORTEDSOURCELUMINOSITYPROBE_HPP

#include "InputModelFormProbe.hpp"
#include "MaterialWavelengthRangeInterface.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedSourceLuminosityProbe probes the specific luminosity aggregated over all imported
    sources in the simulation at the wavelengths specified by a user-configured wavelength grid.
    The probe uses the data as represented by the imported snapshot, without involving the spatial
    grid of the simulation. It produces output only if the simulation has at least one source
    component and if all sources are imported.

    When associated with a form that samples the luminosity at a set of positions, such as for a
    linear or planar cut, the probe outputs a monochromatic luminosity volume density (with SI
    units of W/m/m3 for the per-wavelength flavor). When associated with a form that projects the
    luminosity along a path, the probe outputs a monochromatic luminosity surface density divided
    by the area of the unit sphere (with resulting SI units of W/m/m2/sr for the per-wavelength
    flavor). In case of parallel projection towards a distant observer, this quantity is equivalent
    to surface brightness. In other words, this probe can be used to produce transparent noise-free
    images of the imported sources.

    If the \em convolve configuration property is false, the specific luminosity is sampled at the
    characteristic wavelength of each bin. If the \em convolve property is true, the specific
    luminosity is averaged over the width of each wavelength bin (for disjoint wavelength grids) or
    convolved over each broadband transmission curve (for band wavelength bins). This option is
    ignored for oligochromatic simulation modes. */
class ImportedSourceLuminosityProbe : public InputModelFormProbe, public MaterialWavelengthRangeInterface
{
    ITEM_CONCRETE(ImportedSourceLuminosityProbe, InputModelFormProbe, "imported source: specific luminosity density")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedSourceLuminosityProbe, "ImportedSource")

        PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for this probe")
        ATTRIBUTE_RELEVANT_IF(wavelengthGrid, "Panchromatic")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultInstrumentWavelengthGrid")
        ATTRIBUTE_DISPLAYED_IF(wavelengthGrid, "Level2")

        PROPERTY_BOOL(convolve, "average or convolve the specific luminosity over the wavelength bin")
        ATTRIBUTE_DEFAULT_VALUE(convolve, "false")
        ATTRIBUTE_RELEVANT_IF(convolve, "Panchromatic")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns a wavelength range corresponding to the user-configured wavelength
        grid for this probe, if any, indicating that wavelength-dependent material properties may
        be required for this wavelength range. */
    Range wavelengthRange() const override;

    /** This function returns a pointer to the user-configured wavelength grid for this probe, if
        any, indicating that wavelength-dependent material properties may be required for these
        wavelengths. */
    WavelengthGrid* materialWavelengthGrid() const override;

protected:
    /** This function probes the specified imported source components. */
    void probeImportedSources(const vector<const ImportedSource*>& sources,
                              const vector<const Snapshot*>& snapshots) override;
};

////////////////////////////////////////////////////////////////////

#endif

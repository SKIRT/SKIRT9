/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDMEDIUMTEMPERATUREPROBE_HPP
#define IMPORTEDMEDIUMTEMPERATUREPROBE_HPP

#include "InputModelFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedMediumTemperatureProbe probes the temperature of each imported medium component in the
    simulation for which that information is available. The probe uses the data as represented by
    the imported snapshot, without involving the spatial grid of the simulation. The temperature is
    (mass or number) density-weighted where necessary, i.e. when the probe is associated with a
    form that projects the quantity along a path, or when two or more smoothed particles in the
    imported data overlap.

    For dust medium components, if metallicity is being imported in addition to temperature, the
    weighting is based on the (dusty) gas density rather than the dust density itself. In other
    words, the probe uses the medium density before the automatic multiplication with metallicity
    occurs, but ignoring any cells or particles that do not contain dust (because the corresponding
    information is removed during import). */
class ImportedMediumTemperatureProbe : public InputModelFormProbe
{
    ITEM_CONCRETE(ImportedMediumTemperatureProbe, InputModelFormProbe, "imported medium: temperature")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedMediumTemperatureProbe, "ImportedMedium")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function probes the specified imported medium component. */
    void probeImportedMedium(string sh, const ImportedMedium* medium, const Snapshot* snapshot) override;
};

////////////////////////////////////////////////////////////////////

#endif

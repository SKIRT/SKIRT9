/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDMEDIUMMETALLICITYPROBE_HPP
#define IMPORTEDMEDIUMMETALLICITYPROBE_HPP

#include "InputModelFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedMediumMetallicityProbe probes the metallicity of each imported medium component in the
    simulation for which that information is available. The probe uses the data as represented by
    the imported snapshot, without involving the spatial grid of the simulation. The metallicity is
    (mass or number) density-weighted where necessary, i.e. when the probe is associated with a
    form that projects the quantity along a path, or when two or more smoothed particles in the
    imported data overlap.

    For dust medium components, the weighting is based on the (dusty) gas density rather than the
    dust density itself. In other words, the probe uses the medium density before the automatic
    multiplication with metallicity occurs, but ignoring any cells or particles that do not contain
    dust (because the corresponding information is removed during import). */
class ImportedMediumMetallicityProbe : public InputModelFormProbe
{
    ITEM_CONCRETE(ImportedMediumMetallicityProbe, InputModelFormProbe, "imported medium: metallicity")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedMediumMetallicityProbe, "ImportedMedium")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function probes the specified imported medium component. */
    void probeImportedMedium(string sh, const ImportedMedium* medium, const Snapshot* snapshot) override;
};

////////////////////////////////////////////////////////////////////

#endif

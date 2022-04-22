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
    imported data overlap. For dust medium components, the weighting is based on the gas density
    rather than the dust density, i.e. the density before the automatic multiplication with
    metallicity occurs. */
class ImportedMediumMetallicityProbe : public InputModelFormProbe
{
    ITEM_CONCRETE(ImportedMediumMetallicityProbe, InputModelFormProbe, "the metallicity of the imported medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedMediumMetallicityProbe, "ImportedMedium")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function probes the specified imported medium component. */
    void probeImportedMedium(string sh, const ImportedMedium* medium, const Snapshot* snapshot) override;
};

////////////////////////////////////////////////////////////////////

#endif

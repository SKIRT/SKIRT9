/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDMEDIUMDENSITYPROBE_HPP
#define IMPORTEDMEDIUMDENSITYPROBE_HPP

#include "InputModelFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedMediumDensityProbe probes the density of each imported medium component in the
    simulation as represented by the imported snapshot, without involving the spatial grid of the
    simulation. The probe outputs the mass or number density depending on which one has been
    imported. When associated with a form that projects the quantity along a path, the probe
    outputs mass surface density or number surface density (column density). */
class ImportedMediumDensityProbe : public InputModelFormProbe
{
    ITEM_CONCRETE(ImportedMediumDensityProbe, InputModelFormProbe, "imported medium: mass or number density")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedMediumDensityProbe, "ImportedMedium")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function probes the specified imported medium component. */
    void probeImportedMedium(string sh, const ImportedMedium* medium, const Snapshot* snapshot) override;
};

////////////////////////////////////////////////////////////////////

#endif

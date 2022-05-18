/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDSOURCEMETALLICITYPROBE_HPP
#define IMPORTEDSOURCEMETALLICITYPROBE_HPP

#include "ImportedSourceWeightedProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedSourceMetallicityProbe probes the metallicity of each imported source component in the
    simulation for which that information is available as one of the parameters for the associated
    %SED family. The probe uses the data as represented by the imported snapshot, without involving
    the spatial grid of the simulation. The metallicity is luminosity- or mass-weighted where
    necessary, i.e. when the probe is associated with a form that projects the quantity along a
    path, or when two or more smoothed particles in the imported data overlap. The weighting scheme
    can be user-configured as described for the ImportedSourceWeightedProbe class, from which this
    class derives. */
class ImportedSourceMetallicityProbe : public ImportedSourceWeightedProbe
{
    ITEM_CONCRETE(ImportedSourceMetallicityProbe, ImportedSourceWeightedProbe, "imported source: metallicity")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function probes the imported source component with the specified snapshot and weight
        function. */
    void probeImportedSourceWeighted(string sh, string sweight, const Snapshot* snapshot,
                                     std::function<double(int m)> weight) override;
};

////////////////////////////////////////////////////////////////////

#endif

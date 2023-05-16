/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDSOURCEAGEPROBE_HPP
#define IMPORTEDSOURCEAGEPROBE_HPP

#include "ImportedSourceWeightedProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedSourceAgeProbe probes the average age of the aggregated imported source components in
    the simulation, assuming that this information is available for all sources as one of the
    parameters of the associated %SED family. The probe uses the data as represented by the
    imported snapshot, without involving the spatial grid of the simulation. The age is luminosity-
    or mass-weighted where necessary, i.e. when the probe is associated with a form that projects
    the quantity along a path, or when two or more smoothed particles in the imported data overlap.
    The weighting scheme can be user-configured as described for the ImportedSourceWeightedProbe
    class, from which this class derives.

    The probe produces output only if the simulation has at least one source component, if all
    sources are imported, and if all of these sources offer the age property and the properties
    necessary to perform the requested type of weighting. */
class ImportedSourceAgeProbe : public ImportedSourceWeightedProbe
{
    ITEM_CONCRETE(ImportedSourceAgeProbe, ImportedSourceWeightedProbe, "imported source: age")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function probes the imported source components with the specified snapshots and weight
        function. */
    void probeImportedSourceWeighted(string sweight, const vector<const Snapshot*>& snapshots,
                                     std::function<double(const Snapshot* snapshot, int m)> weight) override;
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INPUTMODELFORMPROBE_HPP
#define INPUTMODELFORMPROBE_HPP

#include "GenericForm.hpp"
#include "Probe.hpp"
class ImportedMedium;
class ImportedSource;
class Snapshot;

////////////////////////////////////////////////////////////////////

/** InputModelFormProbe is a base class for probes that cooperate with any generic Form subclass to
    describe how the considered input model quantity should be probed. This does not include forms
    that require the spatial grid of the simulation to be present. See the ProbeFormBridge class
    for more information.

    This class offers facilities for subclasses that probe imported source or media components.
    Imported sources are always handled together, i.e. aggregated system-wide, because it is not
    meaningful to average properties over each source individually. Imported media are handled per
    component, because aggregating various medium components can be ambiguous (e.g., they may be
    defined in terms of number or mass density). If probes would be developed that access other
    aspects of the input model, this may require additional functionality in this class or an extra
    intermediate layer of subclasses. */
class InputModelFormProbe : public Probe
{
    ITEM_ABSTRACT(InputModelFormProbe, Probe, "an input model form probe")
        ATTRIBUTE_TYPE_DISPLAYED_IF(InputModelFormProbe, "Level2")

        PROPERTY_ITEM(form, GenericForm, "the form describing how this quantity should be probed")
        ATTRIBUTE_DEFAULT_VALUE(form, "ParallelProjectionForm")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is implemented here. It invokes functions defined in the subclass as follows
        (see the class header form the rationale):

        - if \em all source components in the simulation are of type ImportedSource (and there is
        at least one of them), the probeImportedSources() subclass function is called to handle the
        full set of imported sources at the same time.

        - the probeImportedMedium() subclass function is called to handle each medium component of
        type ImportedMedium in the simulation; any other medium components are skipped. */
    void probe() override;

protected:
    /** This function should be implemented by subclasses that probe imported sources. It will be
        called only if all source components in the configured source system are of type
        ImportedSource and there is at least one of them. The source components and corresponding
        snapshots are listed in order of occurrence. The implementation in this base class does
        nothing. */
    virtual void probeImportedSources(const vector<const ImportedSource*>& sources,
                                      const vector<const Snapshot*>& snapshots);

    /** This function should be implemented by subclasses that probe imported media. It will be
        called for each component of type ImportedMedium in the configured source system, in order
        of occurrence. The first argument is a string representation of the zero-based component
        index (including non-imported components). The implementation in this base class does
        nothing. */
    virtual void probeImportedMedium(string sh, const ImportedMedium* medium, const Snapshot* snapshot);
};

////////////////////////////////////////////////////////////////////

#endif

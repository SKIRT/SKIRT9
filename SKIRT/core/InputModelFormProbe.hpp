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
    Specifically, its probe() function calls distinct functions in the subclass for each imported
    source or media component in the simulation. If probes would be developed that access other
    aspects of the input model, then an extra layer of subclasses may be required. */
class InputModelFormProbe : public Probe
{
    ITEM_ABSTRACT(InputModelFormProbe, Probe, "an input model form probe")
        ATTRIBUTE_TYPE_DISPLAYED_IF(InputModelFormProbe, "Level2")

        PROPERTY_ITEM(form, GenericForm, "the form describing how this quantity should be probed")
        ATTRIBUTE_DEFAULT_VALUE(form, "ParallelProjectionForm")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is implemented here. It calls the probeImportedSource() function for each
        component of type ImportedSource in the configured source system, and the
        probeImportedMedium() function for each component of type ImportedMedium in the configured
        medium system. */
    void probe() override;

protected:
    /** This function should be implemented by subclasses that probe imported sources. It will be
        called for each component of type ImportedSource in the configured source system, in order
        of occurrence. The first argument is a string representation of the zero-based component
        index (including non-imported components). The implementation in this base class does
        nothing. */
    virtual void probeImportedSource(string sh, const ImportedSource* source, const Snapshot* snapshot);

    /** This function should be implemented by subclasses that probe imported media. It will be
        called for each component of type ImportedMedium in the configured source system, in order
        of occurrence. The first argument is a string representation of the zero-based component
        index (including non-imported components). The implementation in this base class does
        nothing. */
    virtual void probeImportedMedium(string sh, const ImportedMedium* medium, const Snapshot* snapshot);
};

////////////////////////////////////////////////////////////////////

#endif

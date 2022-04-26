/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDSOURCEDENSITYPROBE_HPP
#define IMPORTEDSOURCEDENSITYPROBE_HPP

#include "InputModelFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedSourceDensityProbe probes the initial or current mass density of each imported source
    component in the simulation for which the information is available. The probe uses the data as
    represented by the imported snapshot, without involving the spatial grid of the simulation.
    When associated with a form that projects the quantity along a path, the probe outputs mass
    surface density.

    In most use cases, one will be interested in probing the current stellar mass of a source
    component. However, this information is usually not automatically available to the SKIRT input
    model. %SED families that model single stellar populations (e.g. BruzualCharlotSEDFamily,
    FSPSSEDFamily) often require the \em initial stellar mass of the population as an input
    parameter, as opposed to the \em current stellar mass. Many other %SED families (e.g.
    MappingsSEDFamily, CastelliKuruczSEDFamily, BlackBodySEDFamily) do not rely (directly) on
    stellar mass information at all.

    Therefore this probe offers the \em massType option. When set to \c InitialMass, the probe will
    use the initial mass imported as one of the %SED family parameters as a suboptimal solution.
    However, assuming that the current mass can be obtained when preparing SKIRT's input data, a
    better solution is to provide it as a separate, additional column to the import process, enable
    the \em importCurrentMass option for the ImportedSource component, and set the \em massType
    option for the probe to \c CurrentMass. Now the probe will use the separately imported current
    mass. In all cases, if a source component does not import the selected mass type, the probe
    produces no output for that component. */
class ImportedSourceDensityProbe : public InputModelFormProbe
{
    /** The enumeration type indicating whether to use initial mass or current mass. */
    ENUM_DEF(MassType, InitialMass, CurrentMass)
        ENUM_VAL(MassType, InitialMass, "initial mass")
        ENUM_VAL(MassType, CurrentMass, "current mass")
    ENUM_END()

    ITEM_CONCRETE(ImportedSourceDensityProbe, InputModelFormProbe, "the mass density of the imported source")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedSourceDensityProbe, "ImportedSource")

        PROPERTY_ENUM(massType, MassType, "type of mass being probed")
        ATTRIBUTE_DEFAULT_VALUE(massType, "InitialMass")
        ATTRIBUTE_DISPLAYED_IF(massType, "CurrentMass")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function probes the specified imported source component. */
    void probeImportedSource(string sh, const ImportedSource* source, const Snapshot* snapshot) override;
};

////////////////////////////////////////////////////////////////////

#endif

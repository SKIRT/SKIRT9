/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDSOURCEWEIGHTEDPROBE_HPP
#define IMPORTEDSOURCEWEIGHTEDPROBE_HPP

#include "InputModelFormProbe.hpp"
#include "MaterialWavelengthRangeInterface.hpp"
#include <functional>

////////////////////////////////////////////////////////////////////

/** ImportedSourceWeightedProbe is an abstract base class for imported source probes that average
    the probed quantity along a path rather than accumulating it. This class offers the options
    that allow the user to configure the weighting mechanism. There are three choices: using the
    specific luminosity at a given wavelength, using initial mass, or using current mass.

    The probes based on this class produce output only if the simulation has at least one source
    component, if all sources are imported, and if all of these sources offer the properties
    necessary to perform the requested type of weighting.

    Luminosity weighting can be accomplished regardless of the configured %SED family as long as
    the source's spectral range includes the specified wavelength.

    For stellar mass weighting, the situation is more complicated. The current stellar mass
    distribution of the source component is usually not automatically available to the SKIRT input
    model. %SED families that model single stellar populations (e.g. BruzualCharlotSEDFamily,
    FSPSSEDFamily) often require the \em initial stellar mass of the population as an input
    parameter, as opposed to the \em current stellar mass. Many other %SED families (e.g.
    MappingsSEDFamily, CastelliKuruczSEDFamily, BlackBodySEDFamily) do not rely (directly) on
    stellar mass information at all.

    Therefore this probe offers two mass weighting options. When set to \c InitialMass, the probe
    will use the initial mass imported as one of the %SED family parameters as a suboptimal
    solution. However, assuming that the current mass can be obtained when preparing SKIRT's input
    data, a better solution is to provide it as a separate, additional column to the import
    process, enable the \em importCurrentMass option for the ImportedSource component, and set the
    \em weight option for the probe to \c CurrentMass. Now the probe will use the separately
    imported current mass. */
class ImportedSourceWeightedProbe : public InputModelFormProbe, public MaterialWavelengthRangeInterface
{
    /** The enumeration type specifying whether to average the probed quantity using luminosity,
        initial mass or current mass. */
    ENUM_DEF(Weight, Luminosity, InitialMass, CurrentMass)
        ENUM_VAL(Weight, Luminosity, "luminosity")
        ENUM_VAL(Weight, InitialMass, "initial mass")
        ENUM_VAL(Weight, CurrentMass, "current mass")
    ENUM_END()

    ITEM_ABSTRACT(ImportedSourceWeightedProbe, InputModelFormProbe, "an imported source weighted probe")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedSourceWeightedProbe, "ImportedSource")

        PROPERTY_ENUM(weight, Weight, "weight for averaging the probed quantity")
        ATTRIBUTE_DEFAULT_VALUE(weight, "Luminosity")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to determine the luminosity")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")
        ATTRIBUTE_RELEVANT_IF(wavelength, "weightLuminosity")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function returns a wavelength range corresponding to the user-configured wavelength,
        indicating that wavelength-dependent material properties may be required for this
        wavelength. */
    Range wavelengthRange() const override;

    //======================== Other Functions =======================

protected:
    /** This function is implemented here. It defines an appropriate weight function and calls the
        probeImportedSourceWeighted() function implemented in the subclass to implement the actual
        probing. */
    void probeImportedSources(const vector<const ImportedSource*>& sources,
                              const vector<const Snapshot*>& snapshots) override;

    /** This function should be implemented by each subclass. The first argument is a string
        representation of the weighting scheme. The last argument is a call-back function that
        returns the weight for the entity with the given index in the snapshot, again according to
        the user configuration. */
    virtual void probeImportedSourceWeighted(string sweight, const vector<const Snapshot*>& snapshots,
                                             std::function<double(const Snapshot* snapshot, int m)> weight) = 0;
};

////////////////////////////////////////////////////////////////////

#endif

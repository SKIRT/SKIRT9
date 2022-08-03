/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDSOURCE_HPP
#define IMPORTEDSOURCE_HPP

#include "Array.hpp"
#include "Range.hpp"
#include "SEDFamily.hpp"
#include "Source.hpp"
class Band;
class Snapshot;

//////////////////////////////////////////////////////////////////////

/** ImportedSource is an abstract class representing a primary radiation source with a spatial and
    spectral luminosity distribution imported from an input file. The input data is usually derived
    from a hydrodynamical simulation snapshot. Various types of snapshots are supported by
    subclasses of this class. Refer to the subclass documentation for information on the file
    format.

    Usually, the input file defines a spatial distribution through smoothed particles, which must
    be interpolated and summed, or through adjacent cells that partition the spatial domain. At the
    level of this abstract class, we use the generic term \em entity for referring to either a
    particle or a cell.

    In addition to spatial information, each entity in the snapshot carries properties that allow
    selecting a particular %SED from a parameterized %SED family. The present class requires the
    user to configure an SEDFamily object for this purpose. The number, type, and order of
    parameters is defined by the %SED family. For each entity, the %SED family is requested to
    select and properly scale a specific %SED based on the entity's properties. Combining the
    spatial and spectral information for an entity yields its contribution to the imported
    radiation source.

    The input file may include a separate column listing the current mass. When this option is
    enabled, the provided current mass can be used for probing the input model. This is relevant
    because %SED families usually do not request the current mass as a parameter (they often use
    the initial mass instead, or do not include direct mass information at all).

    The input file may also include a bulk velocity vector with an optional velocity dispersion for
    each entity. When this option is enabled, the appropriate Doppler shift is taken into account
    when launching photon packets. Apart from the anisotropy resulting from this optional Doppler
    shift, the radiation emitted by this primary source is always isotropic. It is also always
    unpolarized. */
class ImportedSource : public Source
{
    ITEM_ABSTRACT(ImportedSource, Source, "a primary source imported from snapshot data")
        ATTRIBUTE_TYPE_INSERT(ImportedSource, "Dimension3")

        PROPERTY_STRING(filename, "the name of the file to be imported")

        ATTRIBUTE_SUB_PROPERTIES_HERE(ImportedSource)

        PROPERTY_BOOL(importVelocity, "import velocity components (3 columns)")
        ATTRIBUTE_DEFAULT_VALUE(importVelocity, "false")
        ATTRIBUTE_RELEVANT_IF(importVelocity, "Panchromatic")
        ATTRIBUTE_DISPLAYED_IF(importVelocity, "Level2")
        ATTRIBUTE_INSERT(importVelocity, "importVelocity:SourceVelocity")

        PROPERTY_BOOL(importVelocityDispersion, "import velocity dispersion (spherically symmetric)")
        ATTRIBUTE_DEFAULT_VALUE(importVelocityDispersion, "false")
        ATTRIBUTE_RELEVANT_IF(importVelocityDispersion, "Panchromatic&importVelocity")
        ATTRIBUTE_DISPLAYED_IF(importVelocityDispersion, "Level2")

        PROPERTY_BOOL(importCurrentMass, "import current mass")
        ATTRIBUTE_DEFAULT_VALUE(importCurrentMass, "false")
        ATTRIBUTE_DISPLAYED_IF(importCurrentMass, "Level3")
        ATTRIBUTE_INSERT(importCurrentMass, "importCurrentMass:CurrentMass")

        PROPERTY_STRING(useColumns, "a list of names corresponding to columns in the file to be imported")
        ATTRIBUTE_DEFAULT_VALUE(useColumns, "")
        ATTRIBUTE_REQUIRED_IF(useColumns, "false")
        ATTRIBUTE_DISPLAYED_IF(useColumns, "Level3")

        PROPERTY_ITEM(sedFamily, SEDFamily, "the SED family for assigning spectra to the imported sources")
        ATTRIBUTE_DEFAULT_VALUE(sedFamily, "BlackBodySEDFamily")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function imports the snapshot data from the input file through a Snapshot object of
        the appropriate type. Specifically, it first calls the createSnapshot() function, which
        must be implemented in a subclass, to construct and open a Snapshot object of the
        appropriate type. It then passes the user-configurable options of this class to the
        Snapshot object and tells it to import the data.

        Finally, the function constructs a vector with the luminosities (integrated over the
        primary source wavelength range) for all imported entities. This information is used when
        deciding how many photon packets should be launched from each entity. */
    void setupSelfAfter() override;

    /** This function constructs a new Snapshot object of the type appropriate for the subclass,
        calls its open() function, and returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    virtual Snapshot* createAndOpenSnapshot() = 0;

    /** The destructor deletes the snapshot object, if present. */
    ~ImportedSource();

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which is always 3 for an imported
        source. */
    int dimension() const override;

    /** This function returns true if the \em importVelocity flag is enabled for the source. */
    bool hasVelocity() const override;

    /** This function returns the wavelength range for this source. Outside this range, all
        luminosities are zero. This source's wavelength range is determined as the intersection of the
        simulation's source wavelength range (obtained from the simulation configuration) and the
        intrinsic wavelength range of the %SED family associated with the source.

        This function implements the SourceWavelengthRangeInterface interface. */
    Range wavelengthRange() const override;

    /** This function returns the luminosity \f$L\f$ (i.e. radiative power) of the source
        integrated over the wavelength range of primary sources (configured for the source system
        as a whole) and across its complete spatial domain. */
    double luminosity() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
        unit of wavelength) of the source at the specified wavelength \f$\lambda\f$, or zero if the
        wavelength is outside the wavelength range of primary sources (configured for the source
        system as a whole) or if the source simply does not emit at the wavelength. */
    double specificLuminosity(double wavelength) const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ of the source's snapshot
        entity with index \f$m\f$ at the specified wavelength \f$\lambda\f$, or zero if the
        wavelength is outside the wavelength range of primary sources or if the source does not
        emit at the wavelength. If the entity index is out of range, the behavior is undefined.

        This function is intended to provide InputModelProbe instances with access to the
        luminosity per snapshot entity, information that is not otherwise made available to the
        simulation. To preserve proper data encapsulation, this function should \em not be called
        from anywhere else in the simulation machinery. */
    double specificLuminosity(double wavelength, int m) const;

    /** This function returns the average specific luminosity \f$L_\lambda\f$ of the source's
        snapshot entity with index \f$m\f$ in the specified wavelength range, or zero if the
        wavelength range is outside the wavelength range of primary sources or if the source does
        not emit in the wavelength range. If the entity index is out of range, the behavior is
        undefined.

        This function is intended to provide InputModelProbe instances with access to the
        luminosity per snapshot entity, information that is not otherwise made available to the
        simulation. To preserve proper data encapsulation, this function should \em not be called
        from anywhere else in the simulation machinery. */
    double meanSpecificLuminosity(Range wavelengthRange, int m) const;

    /** This function returns the specific luminosity \f$L_\lambda\f$ of the source's snapshot
        entity with index \f$m\f$ convolved over the specified broadband, or zero if the band lies
        outside the wavelength range of primary sources or if the source does not emit in the
        band's wavelength range. If the entity index is out of range, the behavior is undefined.

        This function is intended to provide InputModelProbe instances with access to the
        luminosity per snapshot entity, information that is not otherwise made available to the
        simulation. To preserve proper data encapsulation, this function should \em not be called
        from anywhere else in the simulation machinery. */
    double meanSpecificLuminosity(const Band* band, int m) const;

    /** This function performs some preparations for launching photon packets. It is called in
         serial mode before each segment of photon packet launches, providing the history indices
         mapped by the source system to this particular source. See the description of the
         SourceSystem class for more background information.

         This function distributes the provided range of history indices over the individual
         entities imported by this source, creating a map for use when actually launching the
         photon packets. The number of photon packets allocated to each entity is determined as
         follows:

         \f[ N_m = \left[ (1-\xi) \frac{L_m}{L} + \xi \frac{1}{M} \right] N_s \f]

         where \f$N_s\f$ is the total number of photon packets to be launched by this source,
         \f$N_m\f$ is the number of photon packets to be launched by entity \f$m\f$, \f$L_m\f$ is
         the luminosity of source \f$m\f$, \f$L\f$ is the total luminosity for this source, \f$M\f$
         is the number of entities in this source, and \f$\xi\f$ is the \em emissionBias property
         value of the source system. */
    void prepareForLaunch(double sourceBias, size_t firstIndex, size_t numIndices) override;

    /** This function causes the photon packet \em pp to be launched from the source using the
         given history index and luminosity contribution. It proceeds as follows.

         First, the function finds the entity index that corresponding to the history index using
         the map constructed by the prepareForLaunch() function. It obtains the normalized spectral
         distribition (and the corresponding cumulative distribution) for that entity from the SED
         family configured for this source. In fact, the function sets up a thread-local object
         that caches the spectral distribution for an entity between consecutive invocations of the
         launch() function. This works even if there are multiple sources of this type because each
         thread handles a single photon packet at a time.

         Subsequently, the function samples a wavelength from the entity's SED, properly handling
         the configured wavelength biasing, and asks the Snapshot object to generate a random
         launch position for the entity. If the importVelocity flag is enabled, the function also
         constructs an object that serves a RedshiftInterface appropriate for the velocity of the
         entity. Again, this object is allocated in thread-local storage so that it stays around
         after being handed to the photon packet.

         Finally, the function causes the photon packet to be launched with the information
         described above and an isotropic launch direction. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

    /** This function returns (a pointer to) the snapshot object associated with this imported
        source. It is intended to provide InputModelProbe instances with direct access to the
        snapshot for probing imported information that is not otherwise made available to the
        simulation. To preserve proper data encapsulation, this function should \em not be called
        from anywhere else in the simulation machinery. */
    const Snapshot* snapshot() const;

    //======================== Data Members ========================

private:
    // wavelength information initialized during setup
    bool _oligochromatic{false};                         // true if the simulation is oligochromatic
    Range _wavelengthRange;                              // the wavelength range configured for all primary sources
    double _arbitaryWavelength{0.};                      // an arbitarily chosen wavelength within the source range
    double _xi{0.};                                      // the wavelength bias fraction
    WavelengthDistribution* _biasDistribution{nullptr};  // the wavelength bias distribution

    // snapshot information initialized during setup
    Snapshot* _snapshot{nullptr};
    double _L{0};  // the total bolometric luminosity of all entities (absolute number)
    Array _Lv;     // the relative bolometric luminosity of each entity (normalized to unity)

    // intialized by prepareForLaunch()
    Array _Wv;           // the relative launch weight for each entity (normalized to unity)
    vector<size_t> _Iv;  // first history index allocated to each entity (with extra entry at the end)
};

//////////////////////////////////////////////////////////////////////

#endif

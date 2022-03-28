/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SECONDARYSOURCESYSTEM_HPP
#define SECONDARYSOURCESYSTEM_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
class SecondarySource;
class PhotonPacket;
class ProbePhotonPacketInterface;

//////////////////////////////////////////////////////////////////////

/** SecondarySourceSystem is a helper class that handles much of the complexity involved with
    launching photon packets from secondary sources, i.e. the emission by dust and gas media.

    The class offers an interface that somewhat resembles the interface of the SourceSystem class
    (which holds the primary sources), facilitating the implementation of the photon life cycle for
    both types of sources. Also, bundling this functionality in this class avoids adding everything
    to the already heavy MediumSystem class.

    The SecondarySourceSystem class inherits SimulationItem so that, for example, it can easily use
    the find() template function, but it is not part of the configurable simulation item hierarchy.
    Instead, if secondary emission is enabled in the configuration, the MonteCarloSimulation class
    creates and holds a single instance of the SecondarySourceSystem class. Evidently, the
    SecondarySourceSystem object closely interacts with the MediumSystem object in the
    configuration, which actually holds the media representing the secondary sources, in addition
    to the radiation field obtained by tracing photon packets through these media.

    The functions in this class other than those for constructing and setting up assume that a set
    of photon packets have already been launched for a prior simulation segment, and that radiation
    field information has been accumulated during the life cycles of these packets by calling the
    MediumSystem::storeRadiationField() function. Furthermore, the
    MediumSystem::communicateRadiationField() function must have been called. If this is not the
    case, the behavior of the functions in this class is undefined.

    Media types
    -----------

    The SecondarySourceSystem class supports secondary emission from dust and/or gas medium
    components depending on the configured simulation mode. The \em DustAndGasEmission simulation
    mode allows emission from both media types; the \em DustEmission mode allows emission from dust
    and ignores emission from any other media types; and the \em GasEmission mode allows emission
    from gas and ignores emission from any other media types. In the current implementation,
    electron media types never have any secondary emission.

    For a medium component to support dust emission, it must be configured with a material mix for
    which the MaterialMix::materialType() function returns MaterialType::Dust. The current
    implementation assumes that all material mixes returning this type support secondary emission
    and use the same wavelength grid for discretizing the continuum emission spectrum.
    Specifically, the implementation assumes that for these material mixes, the
    MaterialMix::hasContinuumEmission() function returns true and the
    MaterialMix::emissionWavelengthGrid() function returns the same wavelength grid as the one
    returned by the Configuration::dustEmissionWLG() function. This allows combining all dust
    components into a single secondary source. In other words, for each spatial cell, a single
    aggregate dust emission spectrum is calculated reflecting all dust contents of the cell, and
    photon packets are emitted from this aggregate spectrum.

    For a medium component to support gas emission, it must be configured with a material mix that
    is a subclass of the abstract EmittingGasMix base class and that further adheres to the
    following requirements: the MaterialMix::materialType() function returns MaterialType::Gas (the
    default value implemented by EmittingGasMix) and at least one of the
    MaterialMix::hasContinuumEmission() and MaterialMix::hasLineEmission() functions returns true.
    Each emitting gas component is handled as a separate secondary source; in other words gas
    medium components are not aggregated.

    Spatial cell libraries for dust emission
    ----------------------------------------

    To accelerate processing in cases where the calculation of the secondary emission spectra is
    overly resource intensive, the SecondarySourceSystem class optionally supports a library
    mechanism as described by Baes et al. (2011, ApJS, 196, 22). This library mechanism is limited
    to emission from dust media.

    If so configured by the user, instead of calculating the emission spectrum individually for
    every spatial cell in the system, a library is constructed and template spectra from this
    library are used. Obviously, the templates in the library should be chosen/constructed in such
    a way that they can reasonably approximate the whole range of actual secondary emission spectra
    encountered in the simulation. Different subclasses of the SpatialCellLibrary class achieve
    this goal to different degrees of sophistication. Specifically, a SpatialCellLibrary subclass
    calculates a mapping from each spatial cell \f$m\f$ to the corresponding library entry \f$n\f$
    using a particular heuristic, given the current state of the simulation (i.e. the stored
    radiation field and/or media state).

    \em Important \em note: The calculation of the template spectrum for each library entry assumes
    that all cells mapped to the entry have the same material mix (for each medium component). As a
    result, employing a library scheme other than the AllCellsLibrary (which implements the
    identity mapping and hence essentially no library) does not make much sense in simulations with
    spatially varying material mixes.

    Supporting the library mechanism complicates the procedure described below for distributing
    photon packets over cells. After obtaining the mapping from spatial cells to library entries,
    the prepareForLaunch() function sorts the cells so that all cells mapped to the same library
    entry are consecutive. This allows the launch() function, in turn, to calculate and cache
    information relevant for each library entry and reuse that information for subsequent cells as
    long as they map to the same library entry.

    The launch() function allocates a private DustCellEmission object for each execution thread.
    When presented with a new library entry, this object first determines the average radiation
    field for all spatial cells mapped to the entry. Then, it calls on the DustEmissivity object
    held by the media system to actually calculate the emissivities corresponding to the library
    entry. If the medium system contains multiple dust components \f$h\f$, each with its own
    material mix, the emissivity \f$\varepsilon_{n,h,\ell}\f$ is calculated for each medium
    component \f$h\f$ seperately, and the results are combined into the complete emission spectrum
    for a spatial cell \f$m\f$ through \f[ j_{m,\ell} = \sum_{h=0}^{N_{\text{comp}}-1} \rho_{m,h}
    \, \varepsilon_{n,h,\ell} \f] where \f$\ell\f$ is the wavelength index. Finally, this spectrum
    is normalized to unity. Since the densities \f$\rho_{m,h}\f$ differ for each spatial cell, the
    result must be calculated and stored for each cell separately. If the medium system has only a
    single dust component, the above formula reduces to \f$j_{m,\ell} =\rho_m\,
    \varepsilon_{n,\ell}\f$, so that the normalized emission spectrum is identical for all spatial
    cells that map to a certain library entry.

    Distributing photon packets
    ---------------------------

    One key task of the SecondarySourceSystem object is to distribute photon packet launches across
    the secondary dust and/or gas sources and across the spatial grid used to discretize the media
    in the simulation. In principle, this should/could be achieved by randomly selecting a source
    and a spatial cell for each launch through sampling from an appropriate probability
    distribution. However, similar to the situation with some primary sources (see the
    documentation of the SourceSystem class), a deterministic approach allows significant
    performance optimizations. Because the number of photon packets should be and usually is
    substantially larger than the number of sources and spatial cells, a deterministic approach can
    be considered to be equivalent to the randomized procedure.

    The idea is to iterate through the sources, and for each source through the spatial cells, and
    launch consecutive photon packets from each. The implementation can now construct and cache
    relevant data structures (such as the emission spectrum and the corresponding cumulative
    distribution) for each source/spatial cell combination, and release the information as soon as
    the iteration moves on to the next combination. Because photon packets can (and often are)
    launched in parallel, these data structures must be allocated in thread-local storage, but that
    is only a minor complication.

    For each secondary emission segment (i.e. a sequence of photon packet launches) in the
    simulation, the MonteCarloSimulation object uses the following procedure. It first determines
    the number of photon packets to be launched from the simulation configuration. This number
    \f$N\f$ is passed to the SecondarySourceSystem::prepareForLaunch() function in serial mode.
    Subsequently, the MonteCarloSimulation object launches \f$N\f$ photon packets in (potentially)
    parallel mode, labeling each of the packets with a \em history \em index in the range
    \f$0,...,N-1\f$. While parallel execution threads are working on photon packets in various \em
    chunks of this range, each thread handles photon packets with consecutive history indices
    within a given chunk.

    To achieve the goals described above, the SecondarySourceSystem::prepareForLaunch() function
    maps consecutive history index ranges to each of the source/spatial cell combinations. In a
    first step, the number of photon packets allocated to each secondary source is determined as
    follows:

    \f[ N_s = \left[ (1-\xi) \frac{w_s L_s}{\sum_s w_s L_s} + \xi \frac{w_s}{\sum_s w_s} \right] N
    \f]

    where \f$N\f$ is the total number of photon packets to be launched, \f$N_s\f$ is the number of
    photon packets to be launched by source \f$s\f$, \f$L_s\f$ is the bolometric luminosity of
    source \f$s\f$, \f$w_s\f$ is the \em sourceWeight property value for source \f$s\f$, and
    \f$\xi\f$ is the \em sourceBias property value of the secondary source system. For the purposes
    of this calculation, all dust emission is aggregated into a single source and each emitting gas
    component is handled individually, as explained above. The \em sourceWeight \f$w_s\f$ for the
    aggregated dust emission is configured in DustEmissionOptions, while for each gas component
    this option is configured with the corresponding EmittingGasMix material mix. The \em
    sourceBias \f$\xi\f$ of the secondary source system is configured in SecondaryEmissionOptions.

    Once the photon packets have been distributed across the sources, the number of photon packets
    allocated to each spatial cell for source \f$s\f$ is determined as follows:

    \f[ N_m = \left[ (1-\xi) \frac{L_m}{\sum_m L_m} + \frac{\xi}{M} \right] N_s \f]

    where \f$N_s\f$ is the total number of photon packets to be launched for the source, \f$N_m\f$
    is the number of photon packets to be launched for this source by spatial cell with index
    \f$m\f$, \f$L_m\f$ is the bolometric luminosity for this source of cell \f$m\f$, \f$M\f$ is the
    total number of spatial cells, and \f$\xi\f$ is the \em spatialBias property value configured
    in SecondaryEmissionOptions.
*/
class SecondarySourceSystem : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a SecondarySourceSystem object; it should be invoked only if the
        simulation configuration includes a MediumSystem object and secondary emission is enabled.
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and its
        setup() function has been called. */
    explicit SecondarySourceSystem(SimulationItem* parent);

protected:
    /** This function constructs a secondary source object for each enabled secondary source in the
        simulation, and stores these objects in a list for later use. */
    void setupSelfBefore() override;

public:
    /** This function installs the specified interface as photon packet launch call-back. The
        function probePhotonPacket() provided by the interface will be called for each photon
        packet that is ready to be launched. */
    void installLaunchCallBack(ProbePhotonPacketInterface* callback);

    //======================== Other Functions =======================

public:
    /** This function returns the number of secondary sources in the secondary source system. */
    int numSources() const;

    /** This function prepares the mapping of history indices to sources and tells all sources to
        prepare their individual mapping of history indices to spatial cells. The function returns
        false if the total bolometric luminosity of the secondary sources is zero (which means no
        photon packets can be launched), and true otherwise. */
    bool prepareForLaunch(size_t numPackets);

    /** This function causes the photon packet \em pp to be launched from one of the secondary
        sources, depending on the specified history index. The photon packet's contents is fully
        (re-)initialized so that it is ready to start its lifecycle. */
    void launch(PhotonPacket* pp, size_t historyIndex) const;

    //======================== Data Members ========================

private:
    // initialized by setupSelfBefore()
    vector<SecondarySource*> _sources;  // list of secondary sources managed by this class
    vector<double> _wv;                 // the configured weight for each of these sources
    double _xi;                         // the configured secondary source bias

    // initialized by installLaunchCallBack()
    vector<ProbePhotonPacketInterface*> _callbackv;  // interfaces to be invoked for each packet launch

    // initialized by prepareForLaunch()
    double _L{0};        // the total bolometric luminosity of all sources (absolute number)
    Array _Lv;           // the relative bolometric luminosity of each source (normalized to unity)
    Array _Wv;           // the relative launch weight for each source (normalized to unity)
    double _Lpp{0};      // the average luminosity contribution for each packet
    vector<size_t> _Iv;  // first history index allocated to each source (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif

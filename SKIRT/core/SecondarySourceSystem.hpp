/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SECONDARYSOURCESYSTEM_HPP
#define SECONDARYSOURCESYSTEM_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
class Configuration;
class MediumSystem;
class PhotonPacket;
class ProbePhotonPacketInterface;
class Random;

//////////////////////////////////////////////////////////////////////

/** SecondarySourceSystem is a helper class that handles much of the complexity involved with
    launching photon packets from secondary sources. Important note: in the current version, this
    class supports just emission by dust, not emission by gas.

    The class offers an interface that somewhat resembles the interface of the SourceSystem class
    (which holds the primary sources), facilitating the implementation of the photon life cycle for
    both types of sources. Also, bundling this functionality in this class avoids adding everything
    to the already heavy MediumSystem class.

    The SecondarySourceSystem class inherits SimulationItem, so that it can be located through the
    find() template function, but it is not part of the configurable simulation item hierarchy.
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

    Distributing photon packets
    ---------------------------

    One key task of the SecondarySourceSystem object is to distribute photon packet launches across
    the spatial grid used to discretize the media in the simulation. In principle, this
    should/could be achieved by randomly selecting a spatial cell for each launch through sampling
    from an appropriate probability distribution. However, similar to some primary sources (see the
    documentation of the SourceSystem class), a deterministic approach allows significant
    performance optimizations. Because the number of photon packets should be and usually is
    substantially larger than the number of spatial cells, a deterministic approach can be
    considered to be equivalent to the randomized procedure.

    The idea is to iterate through the spatial cells and launch consecutive photon packets from
    each. The implementation can now construct and cache relevant data structures (such as the
    emission spectrum and the corresponding cumulative distribution) for each spatial cell, and
    release the information as soon as the iteration moves on to the next cell. Because photon
    packets can (and often are) launched in parallel, these data structures must be allocated in
    thread-local storage, but that is only a minor complication.

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
    maps consecutive history index ranges to each of the spatial cells. The number of photon
    packets allocated to each cell is determined as follows:

    \f[ N_m = \left[ (1-\xi) \frac{L_m}{\sum_m L_m} + \frac{\xi}{M} \right] N \f]

    where \f$N\f$ is the total number of photon packets to be launched, \f$N_m\f$ is the number of
    photon packets to be launched by spatial cell with index \f$m\f$, \f$L_m\f$ is the bolometric
    luminosity of cell \f$m\f$, \f$M\f$ is the total number of spatial cells, and \f$\xi\f$ is the
    \em emissionBias property value of the emission simulation mode.

    Spatial cell libraries
    ----------------------

    To accelerate processing in cases where the calculation of the secondary emission spectra is
    overly resource intensive, the SecondarySourceSystem class optionally supports a library
    mechanism as described by Baes et al. (2011, ApJS, 196, 22). If so configured by the user,
    instead of calculating the emission spectrum individually for every spatial cell in the system,
    a library is constructed and template spectra from this library are used. Obviously, the
    templates in the library should be chosen/constructed in such a way that they can reasonably
    approximate the whole range of actual secondary emission spectra encountered in the simulation.
    Different subclasses of the SpatialCellLibrary class achieve this goal to different degrees of
    sophistication. Specifically, a SpatialCellLibrary subclass calculates a mapping from each
    spatial cell \f$m\f$ to the corresponding library entry \f$n\f$ using a particular heuristic,
    given the current state of the simulation (i.e. the stored radiation field and/or media state).

    \em Important \em note: The calculation of the template spectrum for each library entry assumes
    that all cells mapped to the entry have the same material mix (for each medium component). As a
    result, employing a library scheme other than the AllCellsLibrary (which implements the
    identity mapping and hence essentially no library) does not make much sense in simulations with
    spatially varying material mixes.

    Supporting the library mechanism complicates the procedure describe above for distributing
    photon packets. After obtaining the mapping from spatial cells to library entries, the
    prepareForLaunch() function sorts the cells so that all cells mapped to the same library entry
    are consecutive. This allows the launch() function, in turn, to calculate and cache information
    relevant for each library entry and reuse that information for subsequent cells as long as they
    map to the same library entry.

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

    Gas support
    -----------

    To support the continuum emission by gas, the following modifications have been made:

    - prepareForLaunch takes an extra argument, to indicate wether only gas photons, only dust
      photons, or both need to be launched. The ratio is 50/50 for now, but a more customizable
      system might be added later.

    - there is a certain threshold in the list of history indices. Above this threshold, photons
      will be launched using the gas SED instead of the dust SED

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
    /** This function obtains and caches a pointer to several objects in the simulation item
        hierarchy. */
    void setupSelfBefore() override;

public:
    /** This function installs the specified interface as photon packet launch call-back. The
        function probePhotonPacket() provided by the interface will be called for each photon
        packet that is ready to be launched. */
    void installLaunchCallBack(ProbePhotonPacketInterface* callback);

    //======================== Other Functions =======================

public:
    /** This function prepares the mapping of history indices to sources; see the description in
        the class header for more information. The function returns false if the total bolometric
        luminosity of the secondary sources is zero (which means no photon packets can be
        launched), and true otherwise. */
    bool prepareForLaunch(size_t numPackets);

    /** This function causes the photon packet \em pp to be launched from one of the cells in the
        spatial grid using the given history index; see the description in the class header for
        more information. The photon packet's contents is fully (re-)initialized so that it is
        ready to start its lifecycle.

        Before it can randomly emit photon packets from a spatial cell, this function must
        calculate the normalized regular and cumulative dust emission spectrum for the cell from
        the radiation field and the dust properties held by the medium system, using the configured
        emission calculator (LTE or NLTE). Also, it must obtain the average bulk velocity of the
        material in the cell from the medium system. Especially the NLTE emission spectrum
        calculations can be very time-consuming, so it is important to perform them only once for
        each cell. On the other hand, pre-calculating and storing the spectra for all cells in
        advance requires a potentially large amount of memory (proportional to both the number of
        cells and to the number of wavelengths on which the emission is discretized). As described
        in the class header, photon packets launched from a given spatial cell are usually handled
        consecutively by the same execution thread, allowing this function to remember the
        information calculated for the "current" cell from one invocation to the next in a helper
        object allocated with thread-local storage scope. As a result, memory requirements are
        limited to storing the information for only a single cell per execution thread, and the
        calculation is still performed only once per cell.

        Once the emission spectrum for the current cell is known, the function randomly generates a
        wavelength either from this emission spectrum or from the configured bias wavelength
        distribution, adjusting the launch weight with the proper bias factor. It then generates a
        random position uniformly within the spatial cell, and a random direction uniformly on the
        unit sphere (since the emission is assumed to be isotropic). Finally, it actually
        initializes the photon packet with this information. */
    void launch(PhotonPacket* pp, size_t historyIndex) const;

    //======================== Data Members ========================

private:
    // initialized by setupSelfBefore()
    Configuration* _config{nullptr};
    MediumSystem* _ms{nullptr};
    Random* _random{nullptr};

    // initialized by installLaunchCallBack()
    ProbePhotonPacketInterface* _callback{nullptr};  // interface to be invoked for each packet launch if nonzero

    // initialized by prepareForLaunch()
    double _L{0};        // the total bolometric luminosity of all spatial cells
    double _Lpp{0};      // the average luminosity contribution for each packet
    Array _Lv;           // the relative bolometric luminosity of each spatial cell (normalized to unity)
    Array _Wv;           // the relative launch weight for each spatial cell (normalized to unity)
    vector<int> _nv;     // the library entry index corresponding to each spatial cell (i.e. map from cells to entries)
    vector<int> _mv;     // the spatial cell indices sorted so that cells belonging to the same entry are consecutive
    vector<size_t> _Iv;  // first history index allocated to each spatial cell (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif

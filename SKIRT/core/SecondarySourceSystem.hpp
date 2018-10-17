/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SECONDARYSOURCESYSTEM_HPP
#define SECONDARYSOURCESYSTEM_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
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

    The functions in this class other than those for constructing and setting up assume that a set
    of photon packets have been launched for a particular simulation segment, and that radiation
    field information has been accumulated during their life cycles by calling the
    storeRadiationField() function. Furthermore, the communicateRadiationField() function must have
    been called. If this is not the case, the behavior of the functions in this class is undefined.
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
        the class header for more information. The function returns false of the total bolometric
        luminosity of the secondary sources is zero (which means no photon packets can be
        launched), and true otherwise. */
    bool prepareForlaunch(size_t numPackets);

    /** This function causes the photon packet \em pp to be launched from one of the cells in the
        spatial grid using the given history index. The photon packet's contents is fully
        (re-)initialized so that it is ready to start its lifecycle. */
    void launch(PhotonPacket* pp, size_t historyIndex) const;

    //======================== Data Members ========================

private:
    // initialized by setupSelfBefore()
    Configuration* _config{nullptr};
    MediumSystem* _ms{nullptr};
    Random* _random{nullptr};

    // initialized by installLaunchCallBack()
    ProbePhotonPacketInterface* _callback{nullptr}; // interface to be invoked for each packet launch if nonzero

    // initialized by prepareForLaunch()
    double _L{0};       // the total bolometric luminosity of all spatial cells
    double _Lpp{0};     // the average luminosity contribution for each packet
    Array _Lv;          // the relative bolometric luminosity of each spatial cell (normalized to unity)
    Array _Wv;          // the relative launch weight for each spatial cell (normalized to unity)
    vector<size_t> _Iv; // first history index allocated to each spatial cell (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif

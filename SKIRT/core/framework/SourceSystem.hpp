/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SOURCESYSTEM_HPP
#define SOURCESYSTEM_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
#include "Source.hpp"
class PhotonPacket;
class ProbePhotonPacketInterface;

//////////////////////////////////////////////////////////////////////

/** An instance of the SourceSystem class represents a complete primary source system, which is the
    superposition of one or more sources. Each source provides a complete description of its
    radiation, including the spatial and spectral distribution and characteristics such as
    anisotropy and polarization.

    One key task of the SourceSystem object is to distribute photon packet launches across the
    sources. In principle, this should/could be achieved by randomly selecting a source for each
    launch through sampling from an appropriate probability distribution. However, for some
    sources, a deterministic approach allows significant performance optimizations. Because the
    number of photon packets should be and usually is (much) larger than the number of
    (sub)sources, a deterministic approach can be considered to be equivalent to the randomized
    procedure.

    The idea is to iterate through the sources and launch consecutive photon packets from each. A
    source consisting of many subsources (such as particles or cells) can then use a similar
    approach, iterating over these components. The implementation can now construct and cache
    relevant data structures (such as a cumulative spectral distribution) for each subsource, and
    release the information as soon as the iteration moves on to the next subsource. Because photon
    packets can (and often are) launched in parallel, these data structures must be allocated in
    thread-local storage, but that is only a minor complication.

    For each primary emission segment (i.e. a sequence of photon packet launches) in the
    simulation, the MonteCarloSimulation object uses the following procedure. It first determines
    the number of photon packets to be launched from the simulation configuration. This number
    \f$N\f$ is passed to the SourceSystem::prepareForLaunch() function in serial mode.
    Subsequently, the MonteCarloSimulation object launches \f$N\f$ photon packets in (potentially)
    parallel mode, labeling each of the packets with a \em history \em index in the range
    \f$0,...,N-1\f$. While parallel execution threads are working on photon packets in various \em
    chunks of this range, each thread handles photon packets with consecutive history indices
    within a given chunk.

    To achieve the goals described above, the SourceSystem::prepareForLaunch() function maps
    consecutive history index ranges to each of the sources being held. This mapping is also passed
    on to each source, so that it can (but doesn't have to) implement a similar approach for its
    subsources. The number of photon packets allocated to each source is determined as follows:

    \f[ N_s = \left[ (1-\xi) \frac{w_s L_s}{\sum w_s L_s} + \xi \frac{w_s}{\sum w_s} \right] N \f]

    where \f$N\f$ is the total number of photon packets to be launched, \f$N_s\f$ is the number of
    photon packets to be launched by source \f$s\f$, \f$L_s\f$ is the bolometric luminosity of
    source \f$s\f$, \f$w_s\f$ is the \em sourceWeight property value for source \f$s\f$,
    \f$\xi\f$ is the \em sourceBias property value of the source system, and the sums range over
    all sources in the source system.
*/
class SourceSystem : public SimulationItem
{
    ITEM_CONCRETE(SourceSystem, SimulationItem, "a primary source system")

        PROPERTY_DOUBLE(minWavelength, "the shortest wavelength of photon packets launched from primary sources")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(minWavelength, "0.09 micron")
        ATTRIBUTE_RELEVANT_IF(minWavelength, "Panchromatic")

        PROPERTY_DOUBLE(maxWavelength, "the longest wavelength of photon packets launched from primary sources")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(maxWavelength, "100 micron")
        ATTRIBUTE_RELEVANT_IF(maxWavelength, "Panchromatic")

        PROPERTY_DOUBLE_LIST(wavelengths, "the discrete wavelengths of photon packets launched from primary sources")
        ATTRIBUTE_QUANTITY(wavelengths, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelengths, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelengths, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelengths, "0.55 micron")
        ATTRIBUTE_RELEVANT_IF(wavelengths, "Oligochromatic")

        PROPERTY_DOUBLE(sourceBias, "the fraction of photon packets distributed uniformly across primary sources")
        ATTRIBUTE_MIN_VALUE(sourceBias, "[0")
        ATTRIBUTE_MAX_VALUE(sourceBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(sourceBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(sourceBias, "Level3")

        PROPERTY_ITEM_LIST(sources, Source, "the primary sources")
        ATTRIBUTE_DEFAULT_VALUE(sources, "GeometricSource")
        ATTRIBUTE_REQUIRED_IF(sources, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains the bolometric luminosity of each source for later use. */
    void setupSelfAfter() override;

public:
    /** This function installs the specified interface as photon packet launch call-back. The
        function probePhotonPacket() provided by the interface will be called for each photon
        packet that is ready to be launched. */
    void installLaunchCallBack(ProbePhotonPacketInterface* callback);

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source system, which depends on the (lack of)
        symmetry in the geometries of its components. A value of 1 means spherical symmetry, 2
        means axial symmetry and 3 means none of these symmetries. The source with the least
        symmetry (i.e. the highest dimension) determines the result for the whole system. */
    int dimension() const;

    /** This function returns the number of sources in the source system. */
    int numSources() const;

    /** This function returns the bolometric luminosity \f$L\f$ of the source system across its
        spatial and spectral domain, which is the sum of the luminosities of the sources in the
        system. */
    double luminosity() const;

    /** This function prepares the mapping of history indices to sources; see the description in
        the class header for more information. */
    void prepareForLaunch(size_t numPackets);

    /** This function causes the photon packet \em pp to be launched from one of the sources in the
        source system using the given history index. The photon packet's contents is fully
        (re-)initialized so that it is ready to start its lifecycle. */
    void launch(PhotonPacket* pp, size_t historyIndex) const;

    //======================== Data Members ========================

private:
    // intialized during setup
    double _L{0};  // the total bolometric luminosity of all sources (absolute number)
    Array _Lv;     // the relative bolometric luminosity of each source (normalized to unity)
    Array _Wv;     // the relative launch weight for each source (normalized to unity)
    vector<ProbePhotonPacketInterface*> _callbackv;  // interfaces to be invoked for each packet launch

    // intialized by prepareForLaunch()
    double _Lpp{0};      // the average luminosity contribution for each packet
    vector<size_t> _Iv;  // first history index allocated to each source (with extra entry at the end)
};

////////////////////////////////////////////////////////////////

#endif

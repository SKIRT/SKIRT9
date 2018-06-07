/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LAUNCHEDPACKETSPROBE_HPP
#define LAUNCHEDPACKETSPROBE_HPP

#include "Probe.hpp"
#include "ProbePhotonPacketInterface.hpp"
#include "Table.hpp"
#include "WavelengthGrid.hpp"
#include <mutex>

////////////////////////////////////////////////////////////////////

/** LaunchedPacketsProbe outputs a text column file with the number of photon packets launched from
    primary sources on a specified wavelength grid (or on the default instrument wavelength grid).
    The file is named <tt>prefix_launchedpackets.txt</tt>. The first column lists the
    characteristic wavelength of the wavelength bin. The second column lists the total number of
    photon packets launched from primary sources in that wavelength bin (i.e. summed over all
    primary sources). Furthermore, there is an additional column for each primary source, listing
    the number of photon packets launched from that source. The columns are in the same order as
    the sources appear in the configuration file.

    The current implementation uses doubles to count the photon packets in each source/wavelength
    bin. Consequently, the results will be incorrect when the number of photon packets in a single
    bin exceeds 9e15. */
class LaunchedPacketsProbe : public Probe, public ProbePhotonPacketInterface
{
    ITEM_CONCRETE(LaunchedPacketsProbe, Probe, "a probe of the number of photon packets launched from primary sources")

    PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for the photon packet probe")
        ATTRIBUTE_OPTIONAL(wavelengthGrid)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function installs the call-back function for this probe with the source system. */
    void setupSelfBefore() override;

private:
    /** This function resizes and initializes the photon packet counters. It should be called
        exactly once. We cannot do the initialization during setup because Probe subclasses are not
        allowed to cause setup of target objects (such as, in our case, the default instrument
        wavelength grid). */
    void initializeDataMembers();

    //======================== Other Functions =======================

public:
    /** This call-back function increments the photon packet counter for the source that launched
        the given photon packet. */
    void probePhotonPacket(const PhotonPacket* pp) override;

    /** This function outputs the photon packet counts after the simulation run. */
    void probeRun() override;

    //======================== Data Members ========================

private:
    // data members are initialized when the first photon packet is probed
    std::once_flag _initialized;                    // flag to synchronize initialization
    WavelengthGrid* _probeWavelengthGrid{nullptr};  // probe wavelength grid (local or default)
    Table<2> _counts;                               // photon packet counters; indices h, ell
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LAUNCHEDPACKETSPROBE_HPP
#define LAUNCHEDPACKETSPROBE_HPP

#include "ProbePhotonPacketInterface.hpp"
#include "SpecialtyWavelengthGridProbe.hpp"
#include "Table.hpp"

////////////////////////////////////////////////////////////////////

/** LaunchedPacketsProbe outputs a text column file with the number of photon packets launched from
    primary and, if applicable, secondary sources on a specified wavelength grid (or on the default
    instrument wavelength grid). If the simulation iterates over primary and/or secondary emission,
    the photon packets launched during all iterations are accumulated in the counts. The probe uses
    the wavelength at the time when the photon packet was originally emitted, in the rest-frame of
    the original source.

    The output file is named <tt>prefix_launchedpackets.txt</tt>. The first column lists the
    characteristic wavelength of the wavelength bin. Subsequent columns list a number of photon
    packets launched from primary and, if applicable, secondary sources in that wavelength bin, in
    the following order:

    - A column listing the total number of photon packets launched from primary sources.

    - A column for each primary source, listing the number of photon packets launched from that
    source. These columns are in the same order as the source components in the configuration file.

    - If the simulation has secondary emission, a column listing the total number of photon packets
    launched from secondary sources.

    - If the simulation has secondary emission from dust, a column listing the total number of
    photon packets launched from all dust medium components (in other words, all dust emission is
    aggregated).

    - If the simulation has secondary emission from gas, a column for each emitting gas component
    listing the number of photon packets launched from that gas medium component. These columns are
    in the same order as the emitting gas components in the configuration file.

    The current implementation uses doubles to count the photon packets in each source/wavelength
    bin. Consequently, the results will be incorrect when the number of photon packets in a single
    bin exceeds 9e15. */
class LaunchedPacketsProbe : public SpecialtyWavelengthGridProbe, public ProbePhotonPacketInterface
{
    ITEM_CONCRETE(LaunchedPacketsProbe, SpecialtyWavelengthGridProbe,
                  "source: number of photon packets launched from primary and secondary sources")
        ATTRIBUTE_TYPE_DISPLAYED_IF(LaunchedPacketsProbe, "Level2&Source")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns the enumeration \c Run indicating that probing for this probe should
        be performed at the end of the simulation. */
    When when() const override;

    /** This function installs the call-back function for this probe with the primary and secondary
        source systems and initializes the photon packet counters. We do this here as opposed to in
        setupSelfAfter(), because we need the secondary source system (if any) to be constructed
        and fully setup. */
    void initialize() override;

    //======================== Other Functions =======================

public:
    /** This call-back function increments the photon packet counter for the source that launched
        the given photon packet. */
    void probePhotonPacket(const PhotonPacket* pp) override;

protected:
    /** This function outputs the photon packet counts after the simulation run. */
    void probe() override;

    //======================== Data Members ========================

private:
    // data members initialized during initialize
    WavelengthGrid* _probeWavelengthGrid{nullptr};  // probe wavelength grid (local or default)
    Table<2> _primaryCounts;                        // photon packet counters; indices h, ell
    Table<2> _secondaryCounts;                      // photon packet counters; indices h, ell
};

////////////////////////////////////////////////////////////////////

#endif

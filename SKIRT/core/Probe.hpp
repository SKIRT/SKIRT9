/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBE_HPP
#define PROBE_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** Probe is an abstract class representing probes that output information on internal simulation
    data before or during the simulation run. Refer to the ProbeSystem class for more information.

    Probe subclasses \em must adhere to the following rules. In this discussion, \em target refers
    to the item(s) in the simulation hierarchy from which the probe retrieves information.
    - Data encapsulation: a probe can use public interfaces only, even if these interfaces are in
      some cases written specifically to support the probe (refer to the scenarios described in
      the ProbeSystem class).
    - Read-only access: a probe cannot (cause to) change the data structures held by the target.
      There are two exceptions to this rule: (1) during its own setup, a probe can cause setup of
      a target through the find() or interface() functions, and (2) a probe can use a public
      function, provided by the target for this purpose, to install a call-back function that will
      be invoked by the target.
    - Interprobe independence: a probe cannot look for or depend on another probe, nor on the order
      of the various probes in the list held by the probe system.
*/
class Probe : public SimulationItem
{
    ITEM_ABSTRACT(Probe, SimulationItem, "a probe")
    PROPERTY_STRING(probeName, "the name for this probe")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the probe name as human-readable name for the simulation item, so
        that it can be used in log messages to identify the probe and differentiate it from other
        probes. */
    string itemName() const override;

    /** This function is called at the end of the setup phase, i.e. after all simulation items have
        performed setup. It is intended to output relevant information as requested by the user
        configuration. The implementation in this base class does nothing. Each Probe subclass has
        the opportunity to override this function and output something useful. */
    virtual void probeSetup();

    /** This function is called at the end of the run phase, i.e. after all photon packets have
        been emitted and detected. It is intended to output relevant information as requested by
        the user configuration. The implementation in this base class does nothing. Each Probe
        subclass has the opportunity to override this function and output something useful. */
    virtual void probeRun();
};

////////////////////////////////////////////////////////////////////

#endif

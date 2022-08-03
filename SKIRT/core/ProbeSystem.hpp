/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBESYSTEM_HPP
#define PROBESYSTEM_HPP

#include "Probe.hpp"

//////////////////////////////////////////////////////////////////////

/** A ProbeSystem instance maintains a list of zero or more probes (instances of Probe subclasses)
    and helps ensure that these probes are invoked at the appropriate times.

    SKIRT essentially writes two kinds of output. A first type of output can be called "mock
    observations", i.e. the fluxes and surface densities obtained by detecting photon packets
    launched (or peeled-off) during the simulation. This output is written by the Instrument
    subclasses at the end of the simulation run, and it is often the main simulation result. A
    second type of output is information on any of the data structures constructed in preparation
    for or during the simulation run, i.e. internal data. This includes diagnostics used to verify
    configuration and operation of the code and physical quantities that are computed by the
    simulation but cannot be “observed” from the outside. The goal of the probe system is to
    move responsibility for producing the second type of output to a separate set of
    classes/objects as opposed to embedding it in the simulation code itself.

    This design has a number of benefits. From a user perspective, the probe objects can take
    additional attributes to customize the output; for example, requesting a cut through the medium
    density at some offset from the coordinate plane. The probe list can also contain multiple
    probes of the same kind; for example, requesting a cut through the medium density at two or
    more offsets from the coordinate plane. From a developer perspective, there is a more explicit
    interface between the simulation and output code portions, resulting in improved data
    encapsulation. The simulation code is not cluttered with output code, and new probe types can
    be implemented without changing the simulation code (at least, in most cases; see below).

    There are essentially two scenarios for probes to retrieve information from item(s) in the
    simulation hierarchy (the \em target). In the preferred scenario, a probe uses public members
    of the target to retrieve data, and implements the output procedure in the probe. This might
    involve adding some getters to target classes specifically for supporting probe(s). In more
    involved situations, the data can be prepared for output in a public function of the target
    class, while still performing the actual output in the probe, possibly through a callback
    mechanism. This scenario is less preferable because it breaks the clean separation between
    simulation and output code. However, it can be the appropriate choice when the output is best
    implemented through a hierarchy of virtual functions. In any case, probes cannot change the
    data structures held by the target.

    Each Monte Carlo simulation hierarchy has a single ProbeSystem instance, which holds a list
    of Probe subclasses configured by the user. The probeXXX() functions of the probe system,
    and thereby the corresponding probeXXX() functions in each probe object are called at the
    following times during a simulation:

    Name  | Called at the end of
    ------|----------------------
    probeSetup() | the setup phase, i.e. after all simulation items have performed setup
    probeRun()   | the run phase, i.e. after all photon packets have been emitted and detected
    probePrimary()   | each iteration over primary emission
    probeSecondary() | each iteration over secondary emission (which may include primary emission)

    */
class ProbeSystem : public SimulationItem
{
    ITEM_CONCRETE(ProbeSystem, SimulationItem, "a probe system")

        PROPERTY_ITEM_LIST(probes, Probe, "the probes")
        ATTRIBUTE_DEFAULT_VALUE(probes, "!NoMedium:ConvergenceInfoProbe;LuminosityProbe")
        ATTRIBUTE_REQUIRED_IF(probes, "false")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is called at the end of the setup phase, i.e. after all simulation items have
        performed setup. It invokes the function of the same name on all probes in the probe
        system. */
    void probeSetup();

    /** This function is called at the end of the run phase, i.e. after all photon packets have
        been emitted and detected. It invokes the function of the same name on all probes in the
        probe system. */
    void probeRun();

    /** This function is called at the end of each iteration over primary emission, i.e. after all
        photon packets have been processed and the medium state and the radiation field have been
        updated if needed. The function argument specifies the one-based iteration index. This
        function invokes the function of the same name on all probes in the probe system. */
    void probePrimary(int iter);

    /** This function is called at the end of each iteration over secondary emission, i.e. after
        all photon packets have been processed and the medium state and the radiation field have
        been updated if needed. In some execution flows, the iteration may include both a primary
        and secondary emission segment. The function argument specifies the one-based iteration
        index. This function invokes the function of the same name on all probes in the probe
        system. */
    void probeSecondary(int iter);
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBEPHOTONPACKETINTERFACE_HPP
#define PROBEPHOTONPACKETINTERFACE_HPP

#include "Basics.hpp"
class PhotonPacket;

////////////////////////////////////////////////////////////////////

/** ProbePhotonPacketInterface is a pure interface with a single function. Some Probe subclasses
    implement this interface and install it as a call-back with their target object in the
    simulation hierarchy, so that the probePhotonPacket() function gets invoked at predefined times
    in each photon packet's life cycle. */
class ProbePhotonPacketInterface
{
protected:
    /** The empty constructor for the interface. */
    ProbePhotonPacketInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~ProbePhotonPacketInterface() {}

    /** This function probes the given photon packet. After a Probe subclass installs this
        interface as a call-back with its target object, this function gets invoked at predefined
        times in each photon packet's life cycle. Because photon packets may be processed in
        parallel execution threads, this function must be thread-safe. */
    virtual void probePhotonPacket(const PhotonPacket* pp) = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif

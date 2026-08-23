/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EXTINCTIONINTERFACE_HPP
#define EXTINCTIONINTERFACE_HPP

#include "Basics.hpp"
class PhotonPacket;

////////////////////////////////////////////////////////////////////

/** ExtinctionInterface is a pure interface implemented by the MediumSystem class. It offers
    instruments access to the extinction optical depth calculation performed by the medium
    system, without requiring the instrument module to link against the medium system's own
    module. */
class ExtinctionInterface
{
protected:
    /** The empty constructor for the interface. */
    ExtinctionInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~ExtinctionInterface() {}

    /** This function calculates and returns the extinction optical depth along a path through
        the medium system defined by the initial position and direction of the specified
        PhotonPacket object and up to the specified distance. This function is intended for
        handling peel-off photon packets during the photon life cycle. */
    virtual double getExtinctionOpticalDepth(const PhotonPacket* pp, double distance) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

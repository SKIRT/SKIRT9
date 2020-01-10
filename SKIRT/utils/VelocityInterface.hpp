/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VELOCITYINTERFACE_HPP
#define VELOCITYINTERFACE_HPP

#include "Vec.hpp"

////////////////////////////////////////////////////////////////////

/** VelocityInterface is a pure interface to obtain the velocity of a radiation source. */
class VelocityInterface
{
protected:
    /** The empty constructor for the interface. */
    VelocityInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~VelocityInterface() {}

    /** This function returns the velocity of the radiation source. */
    virtual Vec velocity() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

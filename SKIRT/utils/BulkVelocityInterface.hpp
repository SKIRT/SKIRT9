/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BULKVELOCITYINTERFACE_HPP
#define BULKVELOCITYINTERFACE_HPP

#include "Vec.hpp"

////////////////////////////////////////////////////////////////////

/** BulkVelocityInterface is a pure interface to obtain the bulk velocity of a radiation source or
    receiver. */
class BulkVelocityInterface
{
protected:
    /** The empty constructor for the interface. */
    BulkVelocityInterface() { }

public:
    /** The empty destructor for the interface. */
    virtual ~BulkVelocityInterface() { }

    /** This function returns the bulk velocity of the radiation source or receiver. */
    virtual Vec bulkVelocity() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

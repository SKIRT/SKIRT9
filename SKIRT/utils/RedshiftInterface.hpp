/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REDSHIFTINTERFACE_HPP
#define REDSHIFTINTERFACE_HPP

#include "Direction.hpp"

////////////////////////////////////////////////////////////////////

/** RedshiftInterface is a pure interface to obtain the redshift of the radiation emitted by a
    source into a given direction \f$(\theta,\phi)\f$ as a result of the source's bulk velocity. */
class RedshiftInterface
{
protected:
    /** The empty constructor for the interface. */
    RedshiftInterface() { }

public:
    /** The empty destructor for the interface. */
    virtual ~RedshiftInterface() { }

    /** This function returns the redshift of the radiation emitted by a source into the given
        direction \f$(\theta,\phi)\f$ as a result of the source's bulk velocity. For a source that
        is at rest, this function would return zero. */
    virtual double redshiftForDirection(Direction bfk) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

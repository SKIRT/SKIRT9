/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SITELISTINTERFACE_HPP
#define SITELISTINTERFACE_HPP

#include "Position.hpp"

////////////////////////////////////////////////////////////////////

/** SiteListInterface is a pure interface that provides the positions of a set of \em sites. It is
    implemented by density distributions (such as geometries or transfer media) that are defined
    through a set of particles or cells. The locations of these entities can serve, for example, as
    a hint for building an appropriate spatial grid. */
class SiteListInterface
{
protected:
    /** The empty constructor for the interface. */
    SiteListInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~SiteListInterface() {}

    /** This function returns the number of sites in the list. */
    virtual int numSites() const = 0;

    /** This function returns the coordinates of the site with the specified zero-based index. If
        the index is out of range, the behavior is undefined. */
    virtual Position sitePosition(int index) const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif

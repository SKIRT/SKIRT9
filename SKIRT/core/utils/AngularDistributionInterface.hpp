/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ANGULARDISTRIBUTIONINTERFACE_HPP
#define ANGULARDISTRIBUTIONINTERFACE_HPP

#include "Direction.hpp"

////////////////////////////////////////////////////////////////////

/** AngularDistributionInterface is a pure interface to obtain information on an angular
    probability distribution \f$P(\Omega)\,{\mathrm{d}}\Omega\f$. Specifically, the interface
    allows to retrieve the probability corresponding to a given direction \f$(\theta,\phi)\f$. By
    convention, the probability distribution is normalized on the unit sphere as follows: \f[ \int
    P(\Omega) \,{\mathrm{d}}\Omega = \int_{\phi=0}^{2\pi} \int_{\theta=0}^{\pi}
    P(\theta,\phi)\sin\theta \,{\mathrm{d}}\theta \,{\mathrm{d}}\phi = 4\pi\f] */
class AngularDistributionInterface
{
protected:
    /** The empty constructor for the interface. */
    AngularDistributionInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~AngularDistributionInterface() {}

    /** This function returns the probability \f$P(\Omega)\f$ for the given direction
        \f$(\theta,\phi)\f$. For an isotropic distribution, this function would return 1 for any
        direction. */
    virtual double probabilityForDirection(Direction bfk) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

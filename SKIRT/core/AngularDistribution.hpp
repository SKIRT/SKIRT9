/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ANGULARDISTRIBUTION_HPP
#define ANGULARDISTRIBUTION_HPP

#include "AngularDistributionInterface.hpp"
#include "SimulationItem.hpp"
class Random;

////////////////////////////////////////////////////////////////////

/** An instance of an AngularDistribution subclass represents an angular probability distribution
    \f$P(\Omega)\,{\mathrm{d}}\Omega\f$ that serves to describe the anisotropic emission of a point
    source. By convention, the probability distribution is normalized on the unit sphere as
    follows: \f[ \int P(\Omega) \,{\mathrm{d}}\Omega = \int_{\phi=0}^{2\pi} \int_{\theta=0}^{\pi}
    P(\theta,\phi)\sin\theta \,{\mathrm{d}}\theta \,{\mathrm{d}}\phi = 4\pi.\f]

    Subclasses must provide a means to obtain the emission probability \f$P(\Omega)\f$ for a
    given direction \f$(\theta,\phi)\f$, and to generate a random direction \f$(\theta,\phi)\f$
    drawn from the probability distribution \f$P(\Omega)\,{\mathrm{d}}\Omega\f$.

    This class inherits the AngularDistributionInterface offering the probabilityForDirection()
    which much be implemented by a subclass.
*/
class AngularDistribution : public SimulationItem, public AngularDistributionInterface
{
    ITEM_ABSTRACT(AngularDistribution, SimulationItem, "an angular emission profile")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the angular distribution, which depends on its (lack
        of) symmetry. A value of 1 means spherical symmetry, 2 means axial symmetry and 3 means
        none of these symmetries. The function's implementation must be provided in a subclass. */
    virtual int dimension() const = 0;

    /** This function generates a random direction \f$(\theta,\phi)\f$ drawn from the probability
        distribution \f$P(\Omega)\,{\mathrm{d}}\Omega\f$. The function's implementation must be
        provided in a subclass. */
    virtual Direction generateDirection() const = 0;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    //======================== Data Members ========================

private:
    // data member initialized during setup
    Random* _random{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif

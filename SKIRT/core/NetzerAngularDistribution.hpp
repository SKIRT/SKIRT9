/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NETZERANGULARDISTRIBUTION_HPP
#define NETZERANGULARDISTRIBUTION_HPP

#include "AngularDistribution.hpp"
#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The NetzerAngularDistribution class describes anisotropic AGN accretion disk emission as
    proposed by Netzer (1987, MNRAS.225...55N, eq (5)): \f[L(\theta)\propto \begin{cases}
    \cos\theta\,(2\cos\theta+1) & 0\le\theta\le\pi/2 \\ \cos\theta\,(2\cos\theta-1) &
    \pi/2\le\theta\le\pi \end{cases} \f]

    The accretion disk is assumed to be approximated by a point source. The emission pattern is
    axisymmetric, so this angular distribition has a dimension of 2. */
class NetzerAngularDistribution : public AngularDistribution
{
    ITEM_CONCRETE(NetzerAngularDistribution, AngularDistribution, "a Netzer AGN accretion disk emission profile")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function contructs a vector with the cumulative distribution of the anisotropic
        luminosity in function of \f$\theta\f$. For the Netzer luminosity function \f$L(\theta)\f$
        defined in the header of this class, the cumulative distribution is given by \f[ X(\theta)
        \propto \int_0^\theta L(\theta') \sin\theta' \,{\mathrm{d}}\theta' \f] which, after proper
        normalization, leads to \f[ X(\theta) = \begin{cases} \frac{1}{2} - \frac{2}{7}\cos^3\theta
        - \frac{3}{14}\cos^2\theta & 0\le\theta\le\pi/2 \\ \frac{1}{2} - \frac{2}{7}\cos^3\theta +
        \frac{3}{14}\cos^2\theta & \pi/2\le\theta\le\pi \end{cases} \f] */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the angular distribution, which is 2 in this case.
        */
    int dimension() const override;

    /** This function returns the normalized probability for a given direction \f${\bf{k}} =
        (\theta,\phi)\f$ according to the Netzer luminosity profile. It simply implements a
        properly normalized version of the function \f$L(\theta)\f$ defined in the header of this
        class, subject to the normalization \f[ \frac{1}{4\pi} \iint p({\bf{k}})\, {\text{d}}\Omega
        = 1 \f] */
    double probabilityForDirection(Direction bfk) const override;

    /** This function generates a random direction \f${\bf{k}} = (\phi,\theta)\f$ according to the
        Netzer luminosity profile i.e. with \f$\phi\f$ distributed uniformly over the interval
        \f$0\le\phi\le 2\pi\f$, and \f$\theta\f$ determined from the cumulative distribution
        calculated during setup. */
    Direction generateDirection() const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Array _thetav;
    Array _Xv;
};

////////////////////////////////////////////////////////////////////

#endif

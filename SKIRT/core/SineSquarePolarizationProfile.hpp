/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINESQUAREPOLARIZATIONPROFILE_HPP
#define SINESQUAREPOLARIZATIONPROFILE_HPP

#include "PolarizationProfile.hpp"

////////////////////////////////////////////////////////////////////

/** The SineSquarePolarizationProfile class describes a simple axisymmetric polarization profile that
    can be used for testing. The user can configure the direction of the symmetry axis \f$\bf{s}\f$,
    the maximum linear polarization degree \f$P_\mathrm{L}^\mathrm{max}\f$, and the fixed
    polarization angle \f$\gamma\f$. The angle \f$\theta\f$ between the symmetry axis \f$\bf{s}\f$
    and the photon packet's propagation direction \f$\bf{k}\f$ can be easily found from
    \f$\cos\theta=\bf{s}\cdot\bf{k}\f$, assuming both directions are given as unit vectors.
    The Stokes vector for the polarization profile implemented by this class is then defined as:

    \f[\begin{aligned}
    P_\mathrm{L} &= P_\mathrm{L}^\mathrm{max}\,\sin^2 \theta \\
    I &= 1 \\
    Q &= P_\mathrm{L}\,\cos 2\gamma \\
    U &= P_\mathrm{L}\,\sin 2\gamma \\
    V &= 0 \\
    \bf{n} &= \bf{s} \times \bf{k}
    \end{aligned}\f]

    where \f$\bf{n}\f$ is the normal to the reference plane. The reference plane is the plane through
    the symmetry axis and the photon packet propagation direction. This plane is undefined when the
    directions are parallel (i.e. the photon packet is emitted along the positive or negative
    symmetry axis), but in that case the radiation is unpolarized and the reference plane is
    irrelevant. */
class SineSquarePolarizationProfile : public PolarizationProfile
{
    ITEM_CONCRETE(SineSquarePolarizationProfile, PolarizationProfile,
                  "an axysimmetric sine-square polarization emission profile")
        ATTRIBUTE_TYPE_INSERT(SineSquarePolarizationProfile, "Dimension2")

        PROPERTY_DOUBLE(symmetryX, "the direction of the positive symmetry axis, x component")
        ATTRIBUTE_QUANTITY(symmetryX, "length")
        ATTRIBUTE_DEFAULT_VALUE(symmetryX, "0")
        ATTRIBUTE_INSERT(symmetryX, "symmetryX:Dimension3")

        PROPERTY_DOUBLE(symmetryY, "the direction of the positive symmetry axis, y component")
        ATTRIBUTE_QUANTITY(symmetryY, "length")
        ATTRIBUTE_DEFAULT_VALUE(symmetryY, "0")
        ATTRIBUTE_INSERT(symmetryY, "symmetryY:Dimension3")

        PROPERTY_DOUBLE(symmetryZ, "the direction of the positive symmetry axis, z component")
        ATTRIBUTE_QUANTITY(symmetryZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(symmetryZ, "1")

        PROPERTY_DOUBLE(maxPolarizationDegree, "the maximum linear polarization degree")
        ATTRIBUTE_MIN_VALUE(maxPolarizationDegree, "[0")
        ATTRIBUTE_MAX_VALUE(maxPolarizationDegree, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxPolarizationDegree, "0.1")

        PROPERTY_DOUBLE(polarizationAngle, "the linear polarization angle")
        ATTRIBUTE_QUANTITY(polarizationAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(polarizationAngle, "[0 deg")
        ATTRIBUTE_MAX_VALUE(polarizationAngle, "180 deg]")
        ATTRIBUTE_DEFAULT_VALUE(polarizationAngle, "30 deg")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches a normalized version of the direction of the symmetry axis and some
        other values for later use. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the polarization profile, which depends on its (lack
        of) symmetry. The sine-square emission profile is axisymmetric as long as the symmetry axis
    coincides with the positive or negative z-axis. */
    int dimension() const override;

    /** This function returns the Stokes vector defining the polarization state of the radiation
        emitted into the given direction \f$(\theta,\phi)\f$. For the sine-square emission profile,
        this function returns a StokesVector as describe in the header of this class. */
    StokesVector polarizationForDirection(Direction bfk) const override;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Direction _sym;
    double _cos2gamma{0};
    double _sin2gamma{0};
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef AXANGULARDISTRIBUTION_HPP
#define AXANGULARDISTRIBUTION_HPP

#include "AngularDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** AxAngularDistribution is an intermediate abstract class for representing axisymmetric angular
    distributions with an optionally tilted symmetry axis. The class essentially translates the
    directions of the general 3D base class to inclinations for use by subclasses. More
    specifically, it implements all functions required by the AngularDistribution base class, and
    expects subclasses to implement the probabilityForInclination() and generateInclination()
    functions instead. */
class AxAngularDistribution : public AngularDistribution
{
    ITEM_ABSTRACT(AxAngularDistribution, AngularDistribution, "an axisymmetric angular emission profile")
        ATTRIBUTE_TYPE_INSERT(AxAngularDistribution, "Dimension2")

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

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches a normalized version of the direction of the symmetry axis for later
        use. */
    void setupSelfBefore() override;

    //==================== Functions implemented here ===================

public:
    /** This function returns the dimension of the angular distribution, which is 2 if the symmetry
        axis equals the Z-axis, and 3 otherwise. */
    int dimension() const override;

    /** This function returns the normalized probability for a given direction. It determines the
        angle between the symmetry axis and the specified direction, and then calls on the subclass
        to obtain the probability corresponding to this inclination. */
    double probabilityForDirection(Direction bfk) const override;

    /** This function first calls on the subclass to generate a random inclination from its
        distribution, and then generates a 3D direction with that inclination relative to the
        symmetry axis and a random azimuth. */
    Direction generateDirection() const override;

    //=================== Functions to be implemented in subclasses ==================

    /** This function, to be implemented in a subclass, returns the normalized probability for a
        given inclination cosine relative to the symmetry axis of the axisymmetric angular
        distribution. */
    virtual double probabilityForInclinationCosine(double costheta) const = 0;

    /** This function, to be implemented in a subclass, generates a random inclination cosine
        relative to the symmetry axis of the axisymmetric angular distribution. */
    virtual double generateInclinationCosine() const = 0;

    //======================== Data Members ========================

private:
    // data members initialized during setup
    Direction _sym;
};

////////////////////////////////////////////////////////////////////

#endif

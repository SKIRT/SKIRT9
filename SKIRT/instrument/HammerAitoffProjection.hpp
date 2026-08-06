/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HAMMERAITOFFPROJECTION_HPP
#define HAMMERAITOFFPROJECTION_HPP

#include "AllSkyProjection.hpp"

////////////////////////////////////////////////////////////////////

/** The HammerAitoffProjection class implements the Hammer-Aitoff all-sky projection mapping an
    arbitrary direction on the unit sphere to a point in a rectangular frame, and vice versa. */
class HammerAitoffProjection : public AllSkyProjection
{
    ITEM_CONCRETE(HammerAitoffProjection, AllSkyProjection, "the Hammer-Aitoff all-sky projection")
    ITEM_END()

public:
    /** Performs the forward transformation from geographical to rectangular coordinates. */
    virtual void fromGlobeToRectangle(double longitude, double latitude, double& x, double& y) const override;

    /** Performs the backward transformation from rectangular to geographical coordinates. If the
        specified rectangular coordinate values are outside of the possible range for this
        transform, the function returns false and the value of the output variables is undefined.
        */
    virtual bool fromRectangleToGlobe(double x, double y, double& longitude, double& latitude) const override;
};

////////////////////////////////////////////////////////////////////

#endif

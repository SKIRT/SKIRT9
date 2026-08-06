/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ALLSKYPROJECTION_HPP
#define ALLSKYPROJECTION_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** AllSkyProjection is an abstract class representing projections that map an arbitrary direction
    on the unit sphere to a point in a rectangular frame, and vice versa. This base class defines
    an interface for use by client classes, and each subclass implements a specific projection.

    All-sky projections are usually based on geographic map projections designed in the context of
    charting the earth. Different coordinate conventions are used in each field, and unfortunately
    the standard notation confusingly uses the same symbols to mean something different. To avoid
    this confusion, we use the quantity names rather than mathematical symbols in the discussion
    below, and for the argument names in the functions of this class.

    We effectively define three coordinate systems:
      - Rectangular coordinates (x,y) where x and y each range from -1 to 1 (even if the rectangle
        is usually represented two times wider than tall. Some portion near the corners is usually
        unused by the transformation.
      - Spherical coordinates (inclination,azimuth) identical to the coordinates used elsewhere in
        SKIRT, for example as returned by the Direction class. Inclination ranges from 0 to \f$\pi\f$,
        and azimuth from \f$-\pi\f$ to \f$\pi\f$ (relative to the x-axis).
      - Geographical coordinates (longitude, latitude) with Greenwich/equatorial reference.
        Longitude ranges from \f$-\pi\f$ to \f$\pi\f$, and latitude from \f$-\pi/2\f$ to \f$\pi/2\f$

    Including an east-west flip because an all-sky projection looks from inside the sphere, while a
    geographical map looks from outside of the sphere, conversion from spherical to geographical
    coordinates reads:
      -  longitude = -azimuth
      -  latitude = \f$\pi/2\f$ - inclination

    <b>Projection properties</b>

    Projecting a sphere to a rectangle necessarily introduces distortions to distance, angle,
    shape, and/or area. To properly work with the SKIRT instrument flux calibration, an all-sky
    projection must be preserve area ("equal area"), while it may and probably will distort
    distance, angle and shape. Furthermore, in the rectangular coordinates defined above, the
    portion of the rectangle that maps to the sphere (the "used" portion) must be a circle with
    radius one.

    */
class AllSkyProjection : public SimulationItem
{
    ITEM_ABSTRACT(AllSkyProjection, SimulationItem, "an all-sky projection")
    ITEM_END()

    //=================== Convenience Functions =====================

public:
    /** Performs the forward transformation from spherical to rectangular coordinates. This
        function is implemented here in the base class; it calls fromGlobeToRectangle() with
        adjusted coordinates. */
    void fromSphereToRectangle(double inclination, double azimuth, double& x, double& y) const;

    /** Performs the backward transformation from rectangular to spherical coordinates. If the
        specified rectangular coordinate values are outside of the possible range for this
        transform, the function returns false and the value of the output variables is undefined.
        This function is implemented here in the base class; it calls fromRectangleToGlobe() and
        adjusts the resulting coordinates. */
    bool fromRectangleToSphere(double x, double y, double& inclination, double& azimuth) const;

    //=========== Functions to be implemented in subclass ===========

public:
    /** Performs the forward transformation from geographical to rectangular coordinates. The
        implementation must be provided in a subclass. */
    virtual void fromGlobeToRectangle(double longitude, double latitude, double& x, double& y) const = 0;

    /** Performs the backward transformation from rectangular to geographical coordinates. If the
        specified rectangular coordinate values are outside of the possible range for this
        transform, the function returns false and the value of the output variables is undefined.
        The implementation must be provided in a subclass. */
    virtual bool fromRectangleToGlobe(double x, double y, double& longitude, double& latitude) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif

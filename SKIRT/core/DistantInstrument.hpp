/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DISTANTINSTRUMENT_HPP
#define DISTANTINSTRUMENT_HPP

#include "Instrument.hpp"

////////////////////////////////////////////////////////////////////

/** DistantInstrument is an abstract class for representing instruments positioned at a
    sufficiently large distance from the system, so that the observable sky is a plane
    perpendicular to the line of sight rather than a part of a sphere. Consequently
    parallel projection is used and the distance is only important for flux calibration.

    The direction \f${\boldsymbol{k}}_{\text{obs}} = (\theta,\varphi)\f$ towards the instrument is
    specified through the inclination angle \f$\theta\f$ and the azimuth angle \f$\varphi\f$.
    Finally, the instrument can be rotated about the line of sight by specifying its roll angle
    \f$\omega\f$. The table below lists some typical values for these angles, in degrees.

    <TABLE>
    <TR><TD><I>View</I></TD>     <TD>\f$\theta\f$</TD> <TD>\f$\varphi\f$</TD> <TD>\f$\omega\f$</TD> </TR>
    <TR><TD>XY-plane</TD>        <TD>0</TD>   <TD>0</TD>   <TD>90</TD>  </TR>
    <TR><TD>XZ-plane</TD>        <TD>90</TD>  <TD>-90</TD> <TD>0</TD>   </TR>
    <TR><TD>YZ-plane</TD>        <TD>90</TD>  <TD>0</TD>   <TD>0</TD>   </TR>
    <TR><TD>first octant</TD>    <TD>45</TD>  <TD>45</TD>  <TD>0</TD>   </TR>
    </TABLE>
    */
class DistantInstrument : public Instrument
{
    ITEM_ABSTRACT(DistantInstrument, Instrument, "a distant instrument")

    PROPERTY_DOUBLE(distance, "the distance to the system")
        ATTRIBUTE_QUANTITY(distance, "distance")
        ATTRIBUTE_MIN_VALUE(distance, "]0")

    PROPERTY_DOUBLE(inclination, "the inclination angle θ of the detector")
        ATTRIBUTE_QUANTITY(inclination, "posangle")
        ATTRIBUTE_MIN_VALUE(inclination, "0 deg")
        ATTRIBUTE_MAX_VALUE(inclination, "180 deg")
        ATTRIBUTE_DEFAULT_VALUE(inclination, "0 deg")

    PROPERTY_DOUBLE(azimuth, "the azimuth angle φ of the detector")
        ATTRIBUTE_QUANTITY(azimuth, "posangle")
        ATTRIBUTE_MIN_VALUE(azimuth, "-360 deg")
        ATTRIBUTE_MAX_VALUE(azimuth, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(azimuth, "0 deg")

    PROPERTY_DOUBLE(rollAngle, "the roll angle ω of the detector")
        ATTRIBUTE_QUANTITY(rollAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(rollAngle, "-360 deg")
        ATTRIBUTE_MAX_VALUE(rollAngle, "360 deg")
        ATTRIBUTE_DEFAULT_VALUE(rollAngle, "0 deg")
        ATTRIBUTE_DISPLAYED_IF(rollAngle, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function pre-calculates the directions that need to be returned by the instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** Returns the direction towards the observer, expressed in model coordinates. The provided
        photon packet's launching position is not used; it is considered to be very close to the
        coordinate origin from the observer's standpoint, since the distance is sufficiently large.
        */
    Direction bfkobs(const Position& bfr) const override;

    /** Returns the direction along the positive x-axis of the instrument frame, expressed in model
        coordinates. The function applies the inverse instrument transformation to the pixel
        frame's x-axis. The provided photon packet's launching position is not used because the
        orientation of the instrument frame does not depend on it. */
    Direction bfkx(const Position& bfr) const override;

    /** Returns the direction along the positive y-axis of the instrument frame, expressed in model
        coordinates. The function applies the inverse instrument transformation to the pixel
        frame's y-axis. The provided photon packet's launching position is not used because the
        orientation of the instrument frame does not depend on it. */
    Direction bfky(const Position& bfr) const override;

    //======================== Data Members ========================

    // data members derived from the discoverable properties during setup
private:
    Direction _bfkobs;
    Direction _bfkx;
    Direction _bfky;
};

////////////////////////////////////////////////////////////////////

#endif

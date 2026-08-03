/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PERSPECTIVEINSTRUMENT_HPP
#define PERSPECTIVEINSTRUMENT_HPP

#include "HomogeneousTransform.hpp"
#include "Instrument.hpp"

////////////////////////////////////////////////////////////////////

/** The PerspectiveInstrument class implements a full perspective view of the simulated model, with
    arbitrary placement of the viewport outside or inside the model. Only photon packets arriving
    from the front are recorded; light emitted behind the viewport is ignored. For each wavelength
    the instrument outputs the detected surface brightness in every pixel of a given frame as a
    data cube in a FITS file. The instrument does \em not support recording of spatially integrated
    flux densities.

    The perspective instrument is intended mostly for making movies. Each movie frame is generated
    by a separate perspective instrument with the appropriate parameters. */
class PerspectiveInstrument : public Instrument
{
    ITEM_CONCRETE(PerspectiveInstrument, Instrument, "a perspective instrument (mostly for making movies)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PerspectiveInstrument, "Level2")

        PROPERTY_INT(numPixelsX, "the number of viewport pixels in the horizontal direction")
        ATTRIBUTE_MIN_VALUE(numPixelsX, "25")
        ATTRIBUTE_MAX_VALUE(numPixelsX, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsX, "250")

        PROPERTY_INT(numPixelsY, "the number of viewport pixels in the vertical direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "25")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

        PROPERTY_DOUBLE(width, "the width of the viewport")
        ATTRIBUTE_QUANTITY(width, "length")
        ATTRIBUTE_MIN_VALUE(width, "]0")

        PROPERTY_DOUBLE(viewX, "the position of the viewport origin, x component")
        ATTRIBUTE_QUANTITY(viewX, "length")

        PROPERTY_DOUBLE(viewY, "the position of the viewport origin, y component")
        ATTRIBUTE_QUANTITY(viewY, "length")

        PROPERTY_DOUBLE(viewZ, "the position of the viewport origin, z component")
        ATTRIBUTE_QUANTITY(viewZ, "length")

        PROPERTY_DOUBLE(crossX, "the position of the crosshair, x component")
        ATTRIBUTE_QUANTITY(crossX, "length")

        PROPERTY_DOUBLE(crossY, "the position of the crosshair, y component")
        ATTRIBUTE_QUANTITY(crossY, "length")

        PROPERTY_DOUBLE(crossZ, "the position of the crosshair, z component")
        ATTRIBUTE_QUANTITY(crossZ, "length")

        PROPERTY_DOUBLE(upX, "the upwards direction, x component")
        ATTRIBUTE_QUANTITY(upX, "length")

        PROPERTY_DOUBLE(upY, "the upwards direction, y component")
        ATTRIBUTE_QUANTITY(upY, "length")

        PROPERTY_DOUBLE(upZ, "the upwards direction, z component")
        ATTRIBUTE_QUANTITY(upZ, "length")

        PROPERTY_DOUBLE(focal, "the distance from the eye to the viewport origin")
        ATTRIBUTE_QUANTITY(focal, "length")
        ATTRIBUTE_MIN_VALUE(focal, "]0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that all attribute values have been appropriately set and performs
        setup for the instrument. Its most important task is to determine the appropriate
        perspective transformation. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function determines whether the specified instrument has the same observer type,
        position and viewing direction as the receiving instrument, and if so, calls the
        setSameObserverAsPreceding() function to remember the fact. */
    void determineSameObserverAsPreceding(const Instrument* precedingInstrument) override;

    /** Returns the direction towards the eye from the given photon packet launching position. */
    Direction bfkobs(const Position& bfr) const override;

    /** Returns the direction along the positive y-axis of the instrument frame, expressed in model
        coordinates. The provided photon packet's launching position is not used because the
        orientation of the instrument frame does not depend on it. */
    Direction bfky(const Position& bfr) const override;

protected:
    /** This function simulates the detection of a photon packet by the instrument. */
    void detect(PhotonPacket* pp) override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation
    const int& _Nx{_numPixelsX};
    const int& _Ny{_numPixelsY};
    const double& _Sx{_width};
    const double& _Vx{_viewX};
    const double& _Vy{_viewY};
    const double& _Vz{_viewZ};
    const double& _Cx{_crossX};
    const double& _Cy{_crossY};
    const double& _Cz{_crossZ};
    const double& _Ux{_upX};
    const double& _Uy{_upY};
    const double& _Uz{_upZ};
    const double& _Fe{_focal};

    // data members derived from the published attributes during setup
    double _s{0.};                     // width and height of a pixel
    double _Ex{0.}, _Ey{0.}, _Ez{0.};  // eye position
    Direction _bfky;                   // unit vector along the viewport's y-axis
    HomogeneousTransform _transform;   // transform from world to pixel coordinates
};

////////////////////////////////////////////////////////////////////

#endif

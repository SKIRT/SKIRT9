/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HEALPIXSKYINSTRUMENT_HPP
#define HEALPIXSKYINSTRUMENT_HPP

#include "AllSkyProjection.hpp"
#include "HomogeneousTransform.hpp"
#include "Instrument.hpp"

////////////////////////////////////////////////////////////////////

/** Documentation should be written. */
class HEALPixSkyInstrument : public Instrument
{
    ITEM_CONCRETE(HEALPixSkyInstrument, Instrument,
                  "a HEALPix all-sky instrument (for Planck-like observations inside a model)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(AllSkyInstrument, "Level2")

        PROPERTY_ITEM(projection, AllSkyProjection, "the projection used for mapping the sky to a rectangle")
        ATTRIBUTE_DEFAULT_VALUE(projection, "HammerAitoffProjection")

        PROPERTY_INT(numPixelsY, "the number of image pixels in the vertical (shortest) direction")
        ATTRIBUTE_MIN_VALUE(numPixelsY, "25")
        ATTRIBUTE_MAX_VALUE(numPixelsY, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numPixelsY, "250")

        PROPERTY_DOUBLE(radius, "the radius of the observer's all-sky sphere")
        ATTRIBUTE_QUANTITY(radius, "length")
        ATTRIBUTE_MIN_VALUE(radius, "]0")

        PROPERTY_DOUBLE(observerX, "the position of the observer, x component")
        ATTRIBUTE_QUANTITY(observerX, "length")

        PROPERTY_DOUBLE(observerY, "the position of the observer, y component")
        ATTRIBUTE_QUANTITY(observerY, "length")

        PROPERTY_DOUBLE(observerZ, "the position of the observer, z component")
        ATTRIBUTE_QUANTITY(observerZ, "length")

        PROPERTY_DOUBLE(crossX, "the position of the crosshair, x component")
        ATTRIBUTE_QUANTITY(crossX, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossX, "0")

        PROPERTY_DOUBLE(crossY, "the position of the crosshair, y component")
        ATTRIBUTE_QUANTITY(crossY, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossY, "0")

        PROPERTY_DOUBLE(crossZ, "the position of the crosshair, z component")
        ATTRIBUTE_QUANTITY(crossZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(crossZ, "0")

        PROPERTY_DOUBLE(upX, "the upwards direction, x component")
        ATTRIBUTE_QUANTITY(upX, "length")
        ATTRIBUTE_DEFAULT_VALUE(upX, "0")

        PROPERTY_DOUBLE(upY, "the upwards direction, y component")
        ATTRIBUTE_QUANTITY(upY, "length")
        ATTRIBUTE_DEFAULT_VALUE(upY, "0")

        PROPERTY_DOUBLE(upZ, "the upwards direction, z component")
        ATTRIBUTE_QUANTITY(upZ, "length")
        ATTRIBUTE_DEFAULT_VALUE(upZ, "1")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that all attribute values have been appropriately set and performs
        setup for the instrument. Its most important task is to determine an appropriate
        transformation given the instrument position and configuration. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function determines whether the specified instrument has the same observer type,
        position and viewing direction as the receiving instrument, and if so, calls the
        setSameObserverAsPreceding() function to remember the fact. */
    void determineSameObserverAsPreceding(const Instrument* precedingInstrument) override;

    /** Returns the direction towards the observer from the given photon packet launching
        position, expressed in model coordinates. */
    Direction bfkobs(const Position& bfr) const override;

    /** Returns the direction along the positive x-axis in a plane normal to the vector towards the
        observer from the given photon packet launching position, expressed in model coordinates.
        */
    Direction bfkx(const Position& bfr) const override;

    /** Returns the direction along the positive y-axis in a plane normal to the vector towards the
        observer from the given photon packet launching position, expressed in model coordinates.
        */
    Direction bfky(const Position& bfr) const override;

protected:
    /** This function simulates the detection of a photon packet by the instrument. */
    void detect(PhotonPacket* pp) override;

    //======================== Data Members ========================

private:
    // aliases to discoverable data members for ease of notation
    const double& _Ox{_observerX};
    const double& _Oy{_observerY};
    const double& _Oz{_observerZ};
    const double& _Cx{_crossX};
    const double& _Cy{_crossY};
    const double& _Cz{_crossZ};
    const double& _Ux{_upX};
    const double& _Uy{_upY};
    const double& _Uz{_upZ};

    // data members derived from the published attributes during setup
    int _Nx{0};                       // number of pixels in the x direction
    int _Ny{0};                       // number of pixels in the y direction
    double _s{0.};                    // estimated linear size of a pixel
    Direction _bfkx;                  // unit vector along the viewport's x-axis
    Direction _bfky;                  // unit vector along the viewport's y-axis
    HomogeneousTransform _transform;  // transform from world to observer coordinates
};

////////////////////////////////////////////////////////////////////

#endif

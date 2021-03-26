/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SEDINSTRUMENT_HPP
#define SEDINSTRUMENT_HPP

#include "DistantInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** An SEDInstrument object represents a distant instrument that records the spatially integrated
    flux density for each wavelength and outputs an %SED text column file.

    The instrument allows configuring the radius of a circular aperture centered on the origin of
    the model coordinate system and in the plane perpendicular to the instrument's line of sight.
    Photon packets arriving from a point that parallel projects outside of this aperture are
    ignored. If the radius is zero (the default value), the instrument does not have an aperture
    (or, equivalently, the aperture radius is infinite). */
class SEDInstrument : public DistantInstrument
{
    ITEM_CONCRETE(SEDInstrument, DistantInstrument,
                  "a distant instrument that outputs the spatially integrated flux density as an SED")

        PROPERTY_DOUBLE(radius, "the radius of the circular aperture, or zero for no aperture")
        ATTRIBUTE_QUANTITY(radius, "length")
        ATTRIBUTE_MIN_VALUE(radius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(radius, "0")
        ATTRIBUTE_DISPLAYED_IF(radius, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function configures the FluxRecorder instance associated with this instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function simulates the detection of a photon packet by the instrument. It verifies
        that the arriving photon packet projects within the aperture and then calls the detect()
        function of the FluxRecorder instance associated with this instrument. */
    void detect(PhotonPacket* pp) override;

    //======================== Data Members ========================

private:
    // data members derived from the discoverable properties during setup, used in detect()
    double _radius2{0};
    double _costheta{0};
    double _sintheta{0};
    double _cosphi{0};
    double _sinphi{0};
};

////////////////////////////////////////////////////////////////////

#endif

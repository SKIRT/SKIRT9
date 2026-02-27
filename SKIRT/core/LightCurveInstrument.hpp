/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LIGHTCURVEINSTRUMENT_HPP
#define LIGHTCURVEINSTRUMENT_HPP

#include "DistantInstrument.hpp"
#include "TimeGrid.hpp"

////////////////////////////////////////////////////////////////////

/** A LightCurveInstrument object represents a distant instrument that records the spatially
    integrated flux density for each time lag interval and outputs a light curve text column file.

    Although the recorded light curve is spectrally integrated, the instrument still requires a
    wavelength grid to determine its spectral response. Most often, the wavelength grid is used to
    simply limit the spectral range of the recorded photons, but it is also possible to specify an
    arbitrary response curve. Here are a few possibilities:

    - To specify a precise wavelength range with uniform response, use a LinBorderWavelengthGrid
    with a single bin. Note that wavelength grid types without "Border" in their names extrapolate
    borders outside of the specified wavelength range, so it is best to use a "Border" grid.

    - To ensure an identical wavelength range as the one used by other instruments in the
    simulation, simply specify the same wavelength grid. While the multiple bins cause some
    run-time overhead, this is likely to be insignificant.

    - To specify an arbitrary response curve, use a ConfigurableBandWavelengthGrid with one or more
    built-in BroadBand or custom FileBand instances.

    The instrument also allows configuring the radius of a circular aperture centered on the origin
    of the model coordinate system and in the plane perpendicular to the instrument's line of
    sight. Photon packets arriving from a point that parallel projects outside of this aperture are
    ignored. If the radius is zero (the default value), the instrument does not have an aperture
    (or, equivalently, the aperture radius is infinite). */
class LightCurveInstrument : public DistantInstrument
{
    ITEM_CONCRETE(LightCurveInstrument, DistantInstrument,
                  "a distant instrument that outputs the spatially integrated flux density as a light curve")
        ATTRIBUTE_TYPE_ALLOWED_IF(LightCurveInstrument, "ExtinctionOnly")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpectralTimeMapInstrument, "Level2")

        PROPERTY_ITEM(timeGrid, TimeGrid, "the time grid for this instrument")

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

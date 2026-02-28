/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LIGHTCURVEINSTRUMENT_HPP
#define LIGHTCURVEINSTRUMENT_HPP

#include "TimeInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** A LightCurveInstrument object represents a distant instrument with an optional circular
    aperture that records the spatially integrated flux density for each time lag interval and
    outputs a light curve text column file.

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

    */
class LightCurveInstrument : public TimeInstrument
{
    ITEM_CONCRETE(LightCurveInstrument, TimeInstrument,
                  "a distant instrument that outputs the spatially integrated flux density as a light curve")
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
};

////////////////////////////////////////////////////////////////////

#endif

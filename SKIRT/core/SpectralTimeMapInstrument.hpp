/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECTRALTIMEMAPINSTRUMENT_HPP
#define SPECTRALTIMEMAPINSTRUMENT_HPP

#include "ApertureInstrument.hpp"
#include "TimeGrid.hpp"

////////////////////////////////////////////////////////////////////

/** A SpectralTimeMapInstrument object represents a distant instrument with an optional circular
    aperture that records the spatially integrated flux density per wavelength interval and per
    time lag interval, and outputs a FITS file a 2D spectral-time map. */
class SpectralTimeMapInstrument : public ApertureInstrument
{
    ITEM_CONCRETE(SpectralTimeMapInstrument, ApertureInstrument,
                  "a distant instrument that outputs the spatially integrated flux density as a spectral time map")
        ATTRIBUTE_TYPE_ALLOWED_IF(SpectralTimeMapInstrument, "ExtinctionOnly")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpectralTimeMapInstrument, "Level2")

        PROPERTY_ITEM(timeGrid, TimeGrid, "the time grid for this instrument")

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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPECTRALTIMEMAPINSTRUMENT_HPP
#define SPECTRALTIMEMAPINSTRUMENT_HPP

#include "TimeInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** A SpectralTimeMapInstrument object represents a distant instrument with an optional circular
    aperture that records the spatially integrated flux density per wavelength interval and per
    time lag interval, and outputs a FITS file containing a 2D spectral-time map. */
class SpectralTimeMapInstrument : public TimeInstrument
{
    ITEM_CONCRETE(SpectralTimeMapInstrument, TimeInstrument,
                  "a distant instrument that outputs the spatially integrated flux density as a spectral time map")
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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SEDINSTRUMENT_HPP
#define SEDINSTRUMENT_HPP

#include "DistantInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** An SEDInstrument object represents a distant instrument that records the spatially
    integrated flux density for each wavelength and outputs an %SED text column file. */
class SEDInstrument : public DistantInstrument
{
    ITEM_CONCRETE(SEDInstrument, DistantInstrument,
                  "a distant instrument that outputs the spatially integrated flux density as an SED")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function configures the FluxRecorder instance associated with this instrument. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function simulates the detection of a photon packet by the instrument. It calculates
        the optical depth of the path from the photon packet's last interaction site to the
        instrument, and then calls the detect() function of the FluxRecorder instance associated
        with this instrument. */
    void detect(PhotonPacket* pp) override;
};

////////////////////////////////////////////////////////////////////

#endif

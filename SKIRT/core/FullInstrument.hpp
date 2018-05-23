/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FULLINSTRUMENT_HPP
#define FULLINSTRUMENT_HPP

#include "FrameInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** A FullInstrument object represents a distant instrument that records and outputs both the
    spatially integrated flux density for each wavelength (as an %SED text column file) and the
    surface brightness in every pixel of a given frame for each wavelength (as an IFU data cube in
    a FITS file). */
class FullInstrument : public FrameInstrument
{
    ITEM_CONCRETE(FullInstrument, FrameInstrument,
                  "a distant instrument that outputs both the flux density (SED) and surface brightness (data cube)")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function augments the FluxRecorder configuration established in the FrameInstrumet
        base class by requesting an SED in addition to a data cube. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif

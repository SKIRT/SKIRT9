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
    a FITS file).

    This instrument acts as a combination of an SEDInstrument with infinite aperture radius and a
    FrameInstrument with a configurable field of view. Photon packets arriving from a point that
    parallel projects outside of this field of view are ignored for the purposes of the embedded
    FrameInstrument, but are still included in the results of the embedded SEDInstrument.
    Therefore, when interpreting the instrument output for a particular wavelength, spatially
    integrating the surface brightness recorded by the embedded FrameInstrument might produce a
    lower flux density value than the one recorded by the embedded SEDInstrument. */
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

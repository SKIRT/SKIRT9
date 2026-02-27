/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LightCurveInstrument.hpp"
#include "FluxRecorder.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void LightCurveInstrument::setupSelfBefore()
{
    ApertureInstrument::setupSelfBefore();

    // configure flux recorder
    instrumentFluxRecorder()->setTimeGrid(timeGrid());
    instrumentFluxRecorder()->includeLightCurve();
}

////////////////////////////////////////////////////////////////////

void LightCurveInstrument::detect(PhotonPacket* pp)
{
    if (isInsideAperture(pp)) instrumentFluxRecorder()->detect(pp, 0);
}

////////////////////////////////////////////////////////////////////

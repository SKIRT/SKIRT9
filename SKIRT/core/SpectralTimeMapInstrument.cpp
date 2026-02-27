/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpectralTimeMapInstrument.hpp"
#include "FluxRecorder.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void SpectralTimeMapInstrument::setupSelfBefore()
{
    ApertureInstrument::setupSelfBefore();

    // configure flux recorder
    instrumentFluxRecorder()->setTimeGrid(timeGrid());
    instrumentFluxRecorder()->includeSpectralTimeMap();
}

////////////////////////////////////////////////////////////////////

void SpectralTimeMapInstrument::detect(PhotonPacket* pp)
{
    if (isInsideAperture(pp)) instrumentFluxRecorder()->detect(pp, 0);
}

////////////////////////////////////////////////////////////////////

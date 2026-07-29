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
    TimeInstrument::setupSelfBefore();

    // configure flux recorder
    instrumentFluxRecorder()->includeSpectralTimeMap();
}

////////////////////////////////////////////////////////////////////

void SpectralTimeMapInstrument::detect(PhotonPacket* pp)
{
    if (isInsideAperture(pp)) instrumentFluxRecorder()->detect(pp, 0);
}

////////////////////////////////////////////////////////////////////

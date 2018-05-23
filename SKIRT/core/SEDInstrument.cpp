/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SEDInstrument.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void SEDInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    // configure flux recorder
    instrumentFluxRecorder()->includeFluxDensity(distance());
}

////////////////////////////////////////////////////////////////////

void SEDInstrument::detect(PhotonPacket* pp)
{
    //double taupath = opticalDepth(pp);
    double taupath = 0;     // TODO: ask medium system to calculate optical depth
    instrumentFluxRecorder()->detect(pp, 0, taupath);
}

////////////////////////////////////////////////////////////////////

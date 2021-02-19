/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DistantInstrument.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void DistantInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // configure the flux recorder with the appropriate frame and distances
    if (distance() > 0.)
    {
        instrumentFluxRecorder()->setRestFrameDistance(distance());
    }
    else
    {
        auto config = find<Configuration>();
        if (config->redshift() > 0.)
        {
            instrumentFluxRecorder()->setObserverFrameRedshift(config->redshift(), config->angularDiameterDistance(),
                                                               config->luminosityDistance());
        }
        else
        {
            throw FATALERROR("Instrument distance and model redshift are both zero");
        }
    }

    // calculate sine and cosine for our angles
    double costheta = cos(_inclination);
    double sintheta = sin(_inclination);
    double cosphi = cos(_azimuth);
    double sinphi = sin(_azimuth);
    double cosomega = cos(_roll);
    double sinomega = sin(_roll);

    // calculate relevant directions
    _bfkobs = Direction(_inclination, _azimuth);
    _bfky = Direction(-cosphi * costheta * cosomega - sinphi * sinomega,
                      -sinphi * costheta * cosomega + cosphi * sinomega, +sintheta * cosomega);
}

////////////////////////////////////////////////////////////////////

void DistantInstrument::determineSameObserverAsPreceding(const Instrument* precedingInstrument)
{
    auto other = dynamic_cast<const DistantInstrument*>(precedingInstrument);
    if (other && distance() == other->distance() && inclination() == other->inclination()
        && azimuth() == other->azimuth() && roll() == other->roll())
    {
        setSameObserverAsPreceding();
    }
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfkobs(const Position& /*bfr*/) const
{
    return _bfkobs;
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfky(const Position& /*bfr*/) const
{
    return _bfky;
}

////////////////////////////////////////////////////////////////////

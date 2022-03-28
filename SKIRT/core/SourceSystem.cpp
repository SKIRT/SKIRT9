/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SourceSystem.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"
#include "ProbePhotonPacketInterface.hpp"

//////////////////////////////////////////////////////////////////////

void SourceSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // skip preparations if there are no sources
    int Ns = _sources.size();
    if (!Ns) return;

    // obtain the luminosity for each source
    _Lv.resize(Ns);
    for (int h = 0; h != Ns; ++h) _Lv[h] = _sources[h]->luminosity();

    // calculate the total luminosity, and normalize the individual luminosities to unity
    // skip further preparations if the total luminosity is zero
    _L = _Lv.sum();
    if (!_L) return;
    _Lv /= _L;

    // calculate the launch weight for each source, normalized to unity
    Array wv(Ns);
    for (int h = 0; h != Ns; ++h) wv[h] = _sources[h]->sourceWeight();
    Array wLv = wv * _Lv;
    double xi = sourceBias();
    _Wv = (1 - xi) * wLv / wLv.sum() + xi * wv / wv.sum();

    // resize the history index mapping vector
    _Iv.resize(Ns + 1);
}

//////////////////////////////////////////////////////////////////////

void SourceSystem::installLaunchCallBack(ProbePhotonPacketInterface* callback)
{
    _callbackv.push_back(callback);
}

//////////////////////////////////////////////////////////////////////

int SourceSystem::dimension() const
{
    int result = 1;
    for (auto source : _sources) result = max(result, source->dimension());
    return result;
}

//////////////////////////////////////////////////////////////////////

int SourceSystem::numSources() const
{
    return _sources.size();
}

//////////////////////////////////////////////////////////////////////

double SourceSystem::luminosity() const
{
    return _L;
}

//////////////////////////////////////////////////////////////////////

void SourceSystem::prepareForLaunch(size_t numPackets)
{
    if (!_L) throw FATALERROR("Cannot launch primary source photon packets when total luminosity is zero");

    // determine the first history index for each source
    int Ns = _sources.size();
    _Iv[0] = 0;
    double W = 0.;
    for (int h = 1; h != Ns; ++h)
    {
        // track the cumulative normalized weight as a floating point number
        // and limit the index to numPackets to avoid issues with rounding errors
        W += _Wv[h - 1];
        _Iv[h] = min(numPackets, static_cast<size_t>(std::round(W * numPackets)));
    }
    _Iv[Ns] = numPackets;

    //  pass the mapping on to each source
    for (int h = 0; h != Ns; ++h) _sources[h]->prepareForLaunch(sourceBias(), _Iv[h], _Iv[h + 1] - _Iv[h]);

    // calculate the average luminosity contribution for each packet
    _Lpp = _L / numPackets;
}

//////////////////////////////////////////////////////////////////////

void SourceSystem::launch(PhotonPacket* pp, size_t historyIndex) const
{
    // ask the appropriate source to prepare the photon packet for launch
    auto h = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;
    double weight = _Lv[h] / _Wv[h];
    _sources[h]->launch(pp, historyIndex, _Lpp * weight);

    // add additional info
    pp->setPrimaryOrigin(h);

    // invoke any installed launch call-backs
    for (auto callback : _callbackv) callback->probePhotonPacket(pp);
}

//////////////////////////////////////////////////////////////////////

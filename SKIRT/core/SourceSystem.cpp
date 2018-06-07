////////////////////////////////////////////////////////////////////

#include "SourceSystem.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"
#include "ProbePhotonPacketInterface.hpp"

//////////////////////////////////////////////////////////////////////

void SourceSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // obtain the luminosity for each source
    int Ns = _sources.size();
    _Lv.resize(Ns);
    for (int h=0; h!=Ns; ++h) _Lv[h] = _sources[h]->luminosity();

    // calculate the total luminosity, and normalize the individual luminosities to unity
    _L = _Lv.sum();
    _Lv /= _L;

    // calculate the launch weight for each source, normalized to unity
    Array wv(Ns);
    for (int h=0; h!=Ns; ++h) wv[h] = _sources[h]->sourceWeight();
    Array wLv = wv * _Lv;
    double xi = sourceBias();
    _Wv = (1-xi)*wLv/wLv.sum() + xi*wv/wv.sum();

    // resize the history index mapping vector
    _Iv.resize(Ns+1);
}

//////////////////////////////////////////////////////////////////////

void SourceSystem::installLaunchCallBack(ProbePhotonPacketInterface* callback)
{
    if (_callback) throw FATALERROR("Cannot install more than one photon packet launch probe");
    _callback = callback;
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

Range SourceSystem::wavelengthRange() const
{
    return Range(_minWavelength, _maxWavelength);
}

//////////////////////////////////////////////////////////////////////

double SourceSystem::luminosity() const
{
    return _L;
}

//////////////////////////////////////////////////////////////////////

void SourceSystem::prepareForlaunch(size_t numPackets)
{
    // calculate the average luminosity contribution for each packet
    _Lpp = _L / numPackets;

    // determine the first history index for each source and pass the mapping on to each source
    int Ns = _sources.size();
    _Iv[0] = 0;
    for (int h=0; h!=Ns; ++h)
    {
        // limit first index to numPackets to avoid run-over due to rounding errors
        _Iv[h+1] = min(numPackets, _Iv[h] + static_cast<size_t>(std::round(_Wv[h] * numPackets)));
        _sources[h]->prepareForLaunch(_Iv[h], _Iv[h+1]-_Iv[h]);
    }
}

//////////////////////////////////////////////////////////////////////

void SourceSystem::launch(PhotonPacket* pp, size_t historyIndex) const
{
    // ask the appropriate source to prepare the photon packet for launch
    auto h = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;
    double weight = _Lv[h] / _Wv[h];
    _sources[h]->launch(pp, historyIndex, _Lpp*weight);

    // add additional info
    pp->setPrimaryOrigin(h);

    // invoke launch call-back if installed
    if (_callback) _callback->probePhotonPacket(pp);
}

//////////////////////////////////////////////////////////////////////

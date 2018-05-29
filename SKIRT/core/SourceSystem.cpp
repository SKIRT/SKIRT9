////////////////////////////////////////////////////////////////////

#include "SourceSystem.hpp"
#include "NR.hpp"
#include "PhotonPacket.hpp"

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
    for (int h=0; h!=Ns; ++h) wv[h] = _sources[h]->emissionWeight();
    Array wLv = wv * _Lv;
    double xi = emissionBias();
    _Wv = (1-xi)*wLv/wLv.sum() + xi*wv/wv.sum();

    // resize the history index mapping vector
    _Iv.resize(Ns+1);
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
    auto h = std::upper_bound(_Iv.cbegin(), _Iv.cend(), historyIndex) - _Iv.cbegin() - 1;
    double weight = _Lv[h] / _Wv[h];
    _sources[h]->launch(pp, historyIndex, _Lpp*weight);
    pp->setPrimaryOrigin(h);
}

//////////////////////////////////////////////////////////////////////

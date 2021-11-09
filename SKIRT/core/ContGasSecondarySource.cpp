/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ContGasSecondarySource.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

ContGasSecondarySource::ContGasSecondarySource(SimulationItem* parent, int h) : SecondarySource(parent), _h(h) {}

////////////////////////////////////////////////////////////////////

double ContGasSecondarySource::prepareLuminosities()
{
    // avoid warnings
    (void)_h;

    // cache some pointers for later use
    _config = find<Configuration>();
    _ms = find<MediumSystem>();
    _random = find<Random>();

    // calculate and return the total luminosity
    int numCells = _ms->numCells();
    _Lv.resize(numCells);
    _L = _Lv.sum();
    return _L;
}

////////////////////////////////////////////////////////////////////

void ContGasSecondarySource::preparePacketMap(size_t /*firstIndex*/, size_t /*numIndices*/) {}

////////////////////////////////////////////////////////////////////

void ContGasSecondarySource::launch(PhotonPacket* /*pp*/, size_t /*historyIndex*/, double /*L*/) const {}

////////////////////////////////////////////////////////////////////

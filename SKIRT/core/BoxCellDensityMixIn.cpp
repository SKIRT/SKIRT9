/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoxCellDensityMixIn.hpp"
#include "MassInBoxInterface.hpp"
#include "MediumSystem.hpp"

//////////////////////////////////////////////////////////////////////

BoxCellDensityMixIn::BoxCellDensityMixIn() {}

//////////////////////////////////////////////////////////////////////

void BoxCellDensityMixIn::setup(SimulationItem* item)
{
    // verify whether all medium components offer the MassInBoxInterface, without allocating any memory
    auto ms = item->find<MediumSystem>(false);
    if (ms && ms->media().size())
    {
        _enabled = true;
        for (auto medium : ms->media())
            if (!medium->interface<MassInBoxInterface>(0, false)) _enabled = false;
    }

    // if enabled, cache pointer to the interface for each medium component
    if (_enabled)
    {
        for (auto medium : ms->media()) _mibv.push_back(medium->interface<MassInBoxInterface>(0, false));
    }
}

//////////////////////////////////////////////////////////////////////

bool BoxCellDensityMixIn::offersInterface() const
{
    return _enabled;
}

//////////////////////////////////////////////////////////////////////

double BoxCellDensityMixIn::numberDensity(int h, int m) const
{
    Box box = cellBox(m);
    return _mibv[h]->numberInBox(box) / box.volume();
}

//////////////////////////////////////////////////////////////////////

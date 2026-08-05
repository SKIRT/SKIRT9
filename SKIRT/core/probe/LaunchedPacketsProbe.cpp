/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LaunchedPacketsProbe.hpp"
#include "Configuration.hpp"
#include "Indices.hpp"
#include "LockFree.hpp"
#include "PhotonPacket.hpp"
#include "SecondarySourceSystem.hpp"
#include "SourceSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Probe::When LaunchedPacketsProbe::when() const
{
    return When::Run;
}

////////////////////////////////////////////////////////////////////

void LaunchedPacketsProbe::initialize()
{
    // select "local" or default wavelength grid
    _probeWavelengthGrid = find<Configuration>()->wavelengthGrid(wavelengthGrid());
    int numWavelengths = _probeWavelengthGrid->numBins();

    // install ourselves as the launch call-back with the primary source system
    // and resize the primary counts table
    find<SourceSystem>()->installLaunchCallBack(this);
    int numPrimarySources = find<SourceSystem>()->numSources();
    if (numPrimarySources) _primaryCounts.resize(numPrimarySources, numWavelengths);

    // install ourselves as the launch call-back with the secondary source system, if there is one,
    // and resize the secondary counts table, if needed
    auto sss = find<SecondarySourceSystem>(false);
    if (sss)
    {
        sss->installLaunchCallBack(this);
        int numSecondarySources = sss->numSources();
        if (numSecondarySources) _secondaryCounts.resize(numSecondarySources, numWavelengths);
    }
}

////////////////////////////////////////////////////////////////////

void LaunchedPacketsProbe::probePhotonPacket(const PhotonPacket* pp)
{
    // count the packet for each wavelength bin index
    for (int ell : _probeWavelengthGrid->bins(pp->sourceRestFrameWavelength()))
    {
        // get the source component index and register as a primary or secondary packet
        int h = pp->compIndex();
        if (pp->hasPrimaryOrigin())
            LockFree::add(_primaryCounts(h, ell), 1);
        else
            LockFree::add(_secondaryCounts(h, ell), 1);
    }
}

////////////////////////////////////////////////////////////////////

void LaunchedPacketsProbe::probe()
{
    auto units = find<Units>();
    int numWavelengths = _primaryCounts.size(1);
    int numPrimarySources = _primaryCounts.size(0);
    int numSecondarySources = _secondaryCounts.size(0);

    // create a text file and add the columns
    TextOutFile file(this, itemName() + "_launchedpackets", "photon packets launched by primary sources");
    file.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
    file.addColumn("total nr of primary photon packets launched");
    for (int h = 0; h != numPrimarySources; ++h)
        file.addColumn("nr of photon packets launched by primary source " + std::to_string(h + 1));
    if (numSecondarySources)
    {
        file.addColumn("total nr of secondary photon packets launched");
        for (int h = 0; h != numSecondarySources; ++h)
            file.addColumn("nr of photon packets launched by secondary source " + std::to_string(h + 1));
    }

    // write the rows
    for (int ell : Indices(numWavelengths, units->rwavelength()))
    {
        std::vector<double> row;

        // wavelength
        double lambda = _probeWavelengthGrid->wavelength(ell);
        row.push_back(units->owavelength(lambda));

        // primary sources
        {
            size_t primaryTotalIndex = row.size();
            row.push_back(0.);
            for (int h = 0; h != numPrimarySources; ++h)
            {
                row.push_back(_primaryCounts(h, ell));
                row[primaryTotalIndex] += _primaryCounts(h, ell);
            }
        }

        // secondary sources
        if (numSecondarySources)
        {
            size_t secondaryTotalIndex = row.size();
            row.push_back(0.);
            for (int h = 0; h != numSecondarySources; ++h)
            {
                row.push_back(_secondaryCounts(h, ell));
                row[secondaryTotalIndex] += _secondaryCounts(h, ell);
            }
        }

        file.writeRow(row);
    }
}

////////////////////////////////////////////////////////////////////

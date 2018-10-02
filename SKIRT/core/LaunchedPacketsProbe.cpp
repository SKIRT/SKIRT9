/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "LaunchedPacketsProbe.hpp"
#include "Configuration.hpp"
#include "LockFree.hpp"
#include "PhotonPacket.hpp"
#include "SourceSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void LaunchedPacketsProbe::setupSelfAfter()
{
    Probe::setupSelfAfter();

    // install ourselves as the launch call-back with the source system
    find<SourceSystem>()->installLaunchCallBack(this);

    // select "local" or default wavelength grid
    _probeWavelengthGrid = find<Configuration>()->wavelengthGrid(wavelengthGrid());

    // resize the counts table with the appropriate number of sources and number of wavelengths
    int numSources = find<SourceSystem>()->sources().size();
    int numWavelengths = _probeWavelengthGrid->numBins();
    _counts.resize(numSources, numWavelengths);
}

////////////////////////////////////////////////////////////////////

void LaunchedPacketsProbe::probePhotonPacket(const PhotonPacket* pp)
{
    // get the source component index, and abort if this is not a primary source packet
    if (!pp->hasPrimaryOrigin()) return;
    int h = pp->compIndex();

    // count the packet for each wavelength bin index
    for (int ell : _probeWavelengthGrid->bins(pp->sourceRestFrameWavelength()))
        LockFree::add(_counts(h,ell), 1);
}

////////////////////////////////////////////////////////////////////

void LaunchedPacketsProbe::probeRun()
{
    auto units = find<Units>();
    int numSources = _counts.size(0);
    int numWavelengths = _counts.size(1);

    // create a text file and add the columns
    TextOutFile file(this, itemName() + "_launchedpackets", "photon packets launched by primary sources");
    file.addColumn("wavelength", units->uwavelength());
    file.addColumn("total nr of photon packets launched in bin");
    for (int h=0; h!=numSources; ++h)
        file.addColumn("nr of photon packets launched in bin by source " + std::to_string(h+1));

    // write the rows
    for (int ell=0; ell!=numWavelengths; ++ell)
    {
        double lambda = _probeWavelengthGrid->wavelength(ell);

        std::vector<double> row({units->owavelength(lambda), 0.});
        for (int h=0; h!=numSources; ++h)
        {
            row.push_back(_counts(h,ell));
            row[1] += _counts(h,ell);
        }
        file.writeRow(row);
    }
}

////////////////////////////////////////////////////////////////////

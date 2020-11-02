/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileGrainSizeDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileGrainSizeDistribution::setupSelfBefore()
{
    // read the sizes and number of grain values from the input file
    TextInFile infile(this, _filename, "grain size distribution");
    infile.addColumn("grain size", "grainsize", "micron");
    infile.addColumn("size distribution", "pergrainsize", "1/micron");
    infile.readAllColumns(_av, _dndav);
    infile.close();

    // verify the number of values
    if (_av.size() < 2) throw FATALERROR("Number of grain size values is less than 2");
}

////////////////////////////////////////////////////////////////////

double FileGrainSizeDistribution::amin() const
{
    return *begin(_av);
}

////////////////////////////////////////////////////////////////////

double FileGrainSizeDistribution::amax() const
{
    return *(end(_av) - 1);
}

////////////////////////////////////////////////////////////////////

double FileGrainSizeDistribution::dnda(double a) const
{
    return NR::value<NR::interpolateLogLog>(a, _av, _dndav);
}

////////////////////////////////////////////////////////////////////

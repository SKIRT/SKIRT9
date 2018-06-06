/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileWavelengthDistribution.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "Range.hpp"
#include "TextInFile.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

void FileWavelengthDistribution::setupSelfBefore()
{
    WavelengthDistribution::setupSelfBefore();

    // read the wavelengths and probabilities from the input file
    TextInFile infile(this, _filename, "wavelength probability distribution");
    const vector<Array>& columns = infile.readAllColumns(2);
    infile.close();

    // determine the intersected wavelength range
    Range range(columns[0][0]*1e-6, columns[0][columns[0].size()-1]*1e-6);
    range.intersect(interface<WavelengthRangeInterface>()->wavelengthRange());
    if (range.empty()) throw FATALERROR("Wavelength distribution range does not overlap source wavelength range");

    // construct the regular and cumulative distributions in the intersected range
    NR::cdf<NR::interpolateLogLog>(_lambdav, _pv, _Pv, columns[0]*1e-6, columns[1], range);
}

//////////////////////////////////////////////////////////////////////

double FileWavelengthDistribution::probability(double wavelength) const
{
    int i = NR::locateFail(_lambdav, wavelength);
    if (i < 0) return 0.;
    return NR::interpolateLogLog(wavelength, _lambdav[i], _lambdav[i+1], _pv[i], _pv[i+1]);
}

//////////////////////////////////////////////////////////////////////

double FileWavelengthDistribution::generateWavelength() const
{
    return random()->cdf(_lambdav, _Pv);
}

//////////////////////////////////////////////////////////////////////

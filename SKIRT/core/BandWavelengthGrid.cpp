/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BandWavelengthGrid.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void BandWavelengthGrid::setupSelfAfter()
{
    WavelengthGrid::setupSelfAfter();

    // obtain list of bands from subclass
    _bands = bandList();

    // sort the bands in order of pivot wavelength
    std::sort(_bands.begin(), _bands.end(),
              [](Band* b1, Band* b2) { return b1->pivotWavelength() < b2->pivotWavelength(); });

    // verify that no two pivot wavelengths are equal (or very close)
    if (_bands.end() != std::unique(_bands.begin(), _bands.end(), [](Band* b1, Band* b2) {
            return abs(1 - (b1->pivotWavelength() / b2->pivotWavelength())) < 1e-8;
        }))
        throw FATALERROR("Two or more bands have the same pivot wavelength to within 1e-8");
}

////////////////////////////////////////////////////////////////////

int BandWavelengthGrid::numBins() const
{
    return _bands.size();
}

////////////////////////////////////////////////////////////////////

double BandWavelengthGrid::wavelength(int ell) const
{
    return _bands[ell]->pivotWavelength();
}

////////////////////////////////////////////////////////////////////

double BandWavelengthGrid::leftBorder(int ell) const
{
    return _bands[ell]->wavelengthRange().min();
}

////////////////////////////////////////////////////////////////////

double BandWavelengthGrid::rightBorder(int ell) const
{
    return _bands[ell]->wavelengthRange().max();
}

////////////////////////////////////////////////////////////////////

double BandWavelengthGrid::effectiveWidth(int ell) const
{
    return _bands[ell]->effectiveWidth();
}

////////////////////////////////////////////////////////////////////

double BandWavelengthGrid::transmission(int ell, double lambda) const
{
    return _bands[ell]->transmission(lambda);
}

////////////////////////////////////////////////////////////////////

vector<int> BandWavelengthGrid::bins(double lambda) const
{
    vector<int> result;
    int n = _bands.size();
    for (int ell = 0; ell != n; ++ell)
    {
        if (_bands[ell]->wavelengthRange().contains(lambda)) result.push_back(ell);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

int BandWavelengthGrid::bin(double lambda) const
{
    int n = _bands.size();
    for (int ell = 0; ell != n; ++ell)
    {
        if (_bands[ell]->wavelengthRange().contains(lambda)) return ell;
    }
    return -1;
}

////////////////////////////////////////////////////////////////////

const Band* BandWavelengthGrid::band(int ell) const
{
    return _bands[ell];
}

////////////////////////////////////////////////////////////////////

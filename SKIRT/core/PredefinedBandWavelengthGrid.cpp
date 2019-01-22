/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PredefinedBandWavelengthGrid.hpp"
#include "BroadBand.hpp"

////////////////////////////////////////////////////////////////////

vector<Band*> PredefinedBandWavelengthGrid::bandList()
{
    vector<Band*> bands;

    if (_includeGALEX)
    {
        bands.push_back(new BroadBand(this, "GALEX FUV"));
        bands.push_back(new BroadBand(this, "GALEX NUV"));
    }
    if (_includeSDSS)
    {
        bands.push_back(new BroadBand(this, "SDSS u"));
        bands.push_back(new BroadBand(this, "SDSS g"));
        bands.push_back(new BroadBand(this, "SDSS r"));
        bands.push_back(new BroadBand(this, "SDSS i"));
        bands.push_back(new BroadBand(this, "SDSS z"));
    }
    if (_include2MASS)
    {
        bands.push_back(new BroadBand(this, "2MASS J"));
        bands.push_back(new BroadBand(this, "2MASS H"));
        bands.push_back(new BroadBand(this, "2MASS Ks"));
    }
    if (_includeWISE)
    {
        bands.push_back(new BroadBand(this, "WISE W1"));
        bands.push_back(new BroadBand(this, "WISE W2"));
        bands.push_back(new BroadBand(this, "WISE W3"));
        bands.push_back(new BroadBand(this, "WISE W4"));
    }
    if (_includeHERSCHEL)
    {
        bands.push_back(new BroadBand(this, "PACS 70"));
        bands.push_back(new BroadBand(this, "PACS 100"));
        bands.push_back(new BroadBand(this, "PACS 160"));
        bands.push_back(new BroadBand(this, "SPIRE 250"));
        bands.push_back(new BroadBand(this, "SPIRE 350"));
        bands.push_back(new BroadBand(this, "SPIRE 500"));
    }
    return bands;
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedGeometry.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedGeometry::setupSelfAfter()
{
    GenGeometry::setupSelfAfter();

    // create the snapshot with preconfigured mass or density column
    _snapshot = createAndOpenSnapshot();

    // add optional columns if applicable
    if (_importMetallicity) _snapshot->importMetallicity();
    if (_importTemperature) _snapshot->importTemperature();

    // set the density policy;
    // our total mass is normalized to unity, so the constant mass multiplier can have an arbitrary value
    _snapshot->setMassDensityPolicy(1., _importTemperature ? _maxTemperature : 0., true);

    // read the data from file
    _snapshot->readAndClose();

    // remember the normalization factor
    _norm = 1. / _snapshot->mass();
}

//////////////////////////////////////////////////////////////////////

ImportedGeometry::~ImportedGeometry()
{
    delete _snapshot;
}

//////////////////////////////////////////////////////////////////////

double ImportedGeometry::density(Position bfr) const
{
    return _snapshot->density(bfr) * _norm;
}

//////////////////////////////////////////////////////////////////////

Position ImportedGeometry::generatePosition() const
{
    return _snapshot->generatePosition();
}

//////////////////////////////////////////////////////////////////////

double ImportedGeometry::SigmaX() const
{
    return _snapshot->SigmaX() * _norm;
}

//////////////////////////////////////////////////////////////////////

double ImportedGeometry::SigmaY() const
{
    return _snapshot->SigmaY() * _norm;
}

//////////////////////////////////////////////////////////////////////

double ImportedGeometry::SigmaZ() const
{
    return _snapshot->SigmaZ() * _norm;
}

//////////////////////////////////////////////////////////////////////

int ImportedGeometry::numSites() const
{
    return _snapshot->numEntities();
}

//////////////////////////////////////////////////////////////////////

Position ImportedGeometry::sitePosition(int index) const
{
    return _snapshot->position(index);
}

//////////////////////////////////////////////////////////////////////

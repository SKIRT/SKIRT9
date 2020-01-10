/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TabulatedMesh.hpp"
#include "FatalError.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

void TabulatedMesh::setupSelfBefore()
{
    MoveableMesh::setupSelfBefore();

    // get the mesh border points from the subclass
    vector<double> points = getMeshBorderPoints();

    // insert zero point if needed and check basic requirements
    NR::sort(points);
    if (points.size() < 1) throw FATALERROR("The mesh data file has no points");
    if (points.front() < 0.) throw FATALERROR("The mesh data file has negative points");
    if (points.front() != 0.) points.insert(points.begin(), 0.);
    if (points.size() < 2 || points.back() == 0.) throw FATALERROR("The mesh data file has no positive points");

    // assign and scale mesh points
    NR::assign(_mesh, points);
    _mesh /= points.back();

    // set correct number of bins
    setNumBins(_mesh.size() - 1);
}

//////////////////////////////////////////////////////////////////////

Array TabulatedMesh::mesh() const
{
    return _mesh;
}

//////////////////////////////////////////////////////////////////////

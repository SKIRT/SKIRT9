/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileMesh.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "System.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileMesh::setupSelfBefore()
{
    MoveableMesh::setupSelfBefore();

    // read the wavelengths and specific luminosities from the input file
    vector<double> points;
    TextInFile infile(this, _filename, "mesh border points");
    infile.addColumn("border point");
    Array row;
    while (infile.readRow(row)) points.push_back(row[0]);
    infile.close();

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
    setNumBins(_mesh.size()-1);
}

//////////////////////////////////////////////////////////////////////

Array FileMesh::mesh() const
{
    return _mesh;
}

//////////////////////////////////////////////////////////////////////

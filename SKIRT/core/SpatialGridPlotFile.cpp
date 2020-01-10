/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridPlotFile.hpp"
#include "Units.hpp"
#include <iomanip>

////////////////////////////////////////////////////////////////////

SpatialGridPlotFile::SpatialGridPlotFile(const SimulationItem* item, string filename)
    : TextOutFile(item, filename, "data to plot the spatial grid")
{
    // Set the precision
    if (_out.is_open())
    {
        _out << std::setprecision(8);
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeLine(double beg1, double beg2, double end1, double end2)
{
    if (!_out.is_open()) return;

    beg1 = _units->olength(beg1);
    beg2 = _units->olength(beg2);
    end1 = _units->olength(end1);
    end2 = _units->olength(end2);
    _out << beg1 << '\t' << beg2 << '\n';
    _out << end1 << '\t' << end2 << '\n' << '\n';
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeRectangle(double min1, double min2, double max1, double max2)
{
    if (!_out.is_open()) return;

    min1 = _units->olength(min1);
    min2 = _units->olength(min2);
    max1 = _units->olength(max1);
    max2 = _units->olength(max2);
    _out << min1 << '\t' << min2 << '\n';
    _out << min1 << '\t' << max2 << '\n';
    _out << max1 << '\t' << max2 << '\n';
    _out << max1 << '\t' << min2 << '\n';
    _out << min1 << '\t' << min2 << '\n' << '\n';
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeCircle(double radius)
{
    if (!_out.is_open()) return;

    radius = _units->olength(radius);

    for (int l = 0; l <= 360; l++)
    {
        double phi = l * M_PI / 180;
        _out << radius * cos(phi) << '\t' << radius * sin(phi) << '\n';
    }
    _out << '\n';
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeLine(double x1, double y1, double z1, double x2, double y2, double z2)
{
    if (!_out.is_open()) return;

    x1 = _units->olength(x1);
    y1 = _units->olength(y1);
    z1 = _units->olength(z1);
    x2 = _units->olength(x2);
    y2 = _units->olength(y2);
    z2 = _units->olength(z2);
    _out << x1 << '\t' << y1 << '\t' << z1 << '\n' << x2 << '\t' << y2 << '\t' << z2 << '\n' << '\n';
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeCube(double x1, double y1, double z1, double x2, double y2, double z2)
{
    if (!_out.is_open()) return;

    x1 = _units->olength(x1);
    y1 = _units->olength(y1);
    z1 = _units->olength(z1);
    x2 = _units->olength(x2);
    y2 = _units->olength(y2);
    z2 = _units->olength(z2);
    _out << x1 << '\t' << y1 << '\t' << z1 << '\n';
    _out << x2 << '\t' << y1 << '\t' << z1 << '\n';
    _out << x2 << '\t' << y2 << '\t' << z1 << '\n';
    _out << x1 << '\t' << y2 << '\t' << z1 << '\n';
    _out << x1 << '\t' << y1 << '\t' << z1 << '\n' << '\n';

    _out << x1 << '\t' << y1 << '\t' << z2 << '\n';
    _out << x2 << '\t' << y1 << '\t' << z2 << '\n';
    _out << x2 << '\t' << y2 << '\t' << z2 << '\n';
    _out << x1 << '\t' << y2 << '\t' << z2 << '\n';
    _out << x1 << '\t' << y1 << '\t' << z2 << '\n' << '\n';

    _out << x1 << '\t' << y1 << '\t' << z1 << '\n' << x1 << '\t' << y1 << '\t' << z2 << '\n' << '\n';
    _out << x2 << '\t' << y1 << '\t' << z1 << '\n' << x2 << '\t' << y1 << '\t' << z2 << '\n' << '\n';
    _out << x2 << '\t' << y2 << '\t' << z1 << '\n' << x2 << '\t' << y2 << '\t' << z2 << '\n' << '\n';
    _out << x1 << '\t' << y2 << '\t' << z1 << '\n' << x1 << '\t' << y2 << '\t' << z2 << '\n' << '\n';
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writePolyhedron(const vector<double>& coords, const vector<int>& indices)
{
    if (!_out.is_open()) return;

    unsigned int k = 0;
    while (k < indices.size())
    {
        int numvertices = indices[k++];
        int firstindex = indices[k];
        for (int i = 0; i < numvertices; i++)
        {
            int currentindex = indices[k++];
            double x = _units->olength(coords[3 * currentindex]);
            double y = _units->olength(coords[3 * currentindex + 1]);
            double z = _units->olength(coords[3 * currentindex + 2]);
            _out << x << '\t' << y << '\t' << z << '\n';
        }
        double x = _units->olength(coords[3 * firstindex]);
        double y = _units->olength(coords[3 * firstindex + 1]);
        double z = _units->olength(coords[3 * firstindex + 2]);
        _out << x << '\t' << y << '\t' << z << '\n' << '\n';
    }
}

////////////////////////////////////////////////////////////////////

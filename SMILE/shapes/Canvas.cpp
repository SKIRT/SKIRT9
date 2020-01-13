/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Canvas.hpp"
#include "ConvexPolygon.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

Canvas::Canvas(size_t numPixels) : _pixels(numPixels, numPixels, 3), _n(numPixels) {}

////////////////////////////////////////////////////////////////////

void Canvas::setColor(double r, double g, double b)
{
    _r = r;
    _g = g;
    _b = b;
}

////////////////////////////////////////////////////////////////////

void Canvas::setWidth(double w)
{
    _w2 = 0.5 * w;
}

////////////////////////////////////////////////////////////////////

void Canvas::fillBox(double x1, double y1, double x2, double y2)
{
    size_t i1, i2, j1, j2;
    if (indicesInRange(x1, x2, i1, i2) && indicesInRange(y1, y2, j1, j2))
    {
        for (size_t i = i1; i <= i2; ++i)
        {
            for (size_t j = j1; j <= j2; ++j)
            {
                _pixels(i, j, 0) = _r;
                _pixels(i, j, 1) = _g;
                _pixels(i, j, 2) = _b;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void Canvas::fillConvexPolygon(const ConvexPolygon& poly)
{
    size_t j1, j2;
    if (indicesInRange(poly.bottom(), poly.top(), j1, j2))
    {
        for (size_t j = j1; j <= j2; ++j)
        {
            double y = (j + 0.5) / _n;
            size_t i1, i2;
            if (indicesInRange(poly.leftFor(y), poly.rightFor(y), i1, i2))
            {
                for (size_t i = i1; i <= i2; ++i)
                {
                    _pixels(i, j, 0) = _r;
                    _pixels(i, j, 1) = _g;
                    _pixels(i, j, 2) = _b;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void Canvas::drawLine(double x1, double y1, double x2, double y2)
{
    // handle horizontal and vertical cases seperately to avoid singularities
    constexpr double eps = 1e-9;  // this is about 1/n**2 for largest supported n
    double dx = x2 - x1;
    double dy = y2 - y1;
    if (abs(dx) < eps || abs(dy) < eps)
        return fillBox(min(x1, x2) - _w2, min(y1, y2) - _w2, max(x1, x2) + _w2, max(y1, y2) + _w2);

    // determine the projections on the coordinate axes of w/2 along the line
    double l = std::sqrt(dx * dx + dy * dy);
    double wx = _w2 * dx / l;
    double wy = _w2 * dy / l;

    // determine the corner points of the slanted rectangle defining the line boundaries
    // and add them to a polygon object in clockwise order
    ConvexPolygon poly;
    poly.add(x1 - wx + wy, y1 - wx - wy);
    poly.add(x1 - wx - wy, y1 + wx - wy);
    poly.add(x2 + wx - wy, y2 + wx + wy);
    poly.add(x2 + wx + wy, y2 - wx + wy);

    // fill the line boundaries
    fillConvexPolygon(poly);
}

////////////////////////////////////////////////////////////////////

bool Canvas::saveToTiff(string filepath) const
{
    // open the output file
    auto file = System::ofstream(filepath);
    if (file)
    {
        // construct tiff header as uint16 words
        uint16_t n = static_cast<uint16_t>(_n);
        uint16_t lsBytes = static_cast<uint16_t>(_pixels.size() & 65535);
        uint16_t msBytes = static_cast<uint16_t>(_pixels.size() >> 16);
        std::vector<uint16_t> header{{0x4949, 42, 8, 0,      //   0: TIFF header (little endian)
                                      12,                    //   8: number of directory entries
                                                             //  (directory entry: tag,type,count,0,value/offset x 2)
                                      256, 4, 1, 0, n, 0,    //  10: ImageWidth, 1 LONG
                                      257, 4, 1, 0, n, 0,    //  22: ImageLength, 1 LONG
                                      258, 3, 3, 0, 158, 0,  //  34: BitsPerSample, 3 SHORT (-> offset!)
                                      259, 3, 1, 0, 1, 0,    //  46: Compression, 1 SHORT
                                      262, 3, 1, 0, 2, 0,    //  58: PhotometricInterpretation, 1 SHORT
                                      273, 4, 1, 0, 180, 0,  //  70: StripOffsets, 1 LONG
                                      277, 3, 1, 0, 3, 0,    //  82: SamplesPerPixel, 1 SHORT
                                      278, 4, 1, 0, n, 0,    //  94: RowsPerStrip, 1 LONG
                                      279, 4, 1, 0, lsBytes, msBytes,  // 106: StripByteCounts, 1 LONG
                                      282, 5, 1, 0, 164, 0,            // 118: XResolution, 1 RATIONAL (-> offset!)
                                      283, 5, 1, 0, 172, 0,            // 130: YResolution, 1 RATIONAL (-> offset!)
                                      296, 3, 1, 0, 2, 0,              // 142: ResolutionUnit, 1 SHORT
                                      0, 0,                            // 154: IFD list terminator
                                      8, 8, 8,                         // 158: BitsPerSample value
                                      72, 0, 1, 0,                     // 164: XResolution value
                                      72, 0, 1, 0}};                   // 172: YResolution value
                                                                       // 180: Image data
        // write tiff header as bytes
        file.write(reinterpret_cast<char*>(&header[0]), 2 * header.size());

        // write rgb data, 1 byte per color component
        for (size_t jpos = 0; jpos != _n; ++jpos)
        {
            size_t j = _n - jpos - 1;
            for (size_t i = 0; i != _n; ++i)
            {
                for (size_t c = 0; c != 3; ++c)
                {
                    long color = static_cast<long>(_pixels(i, j, c) * 256.);
                    if (color < 0) color = 0;
                    if (color > 255) color = 255;
                    file.put(static_cast<char>(color));
                }
            }
        }
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

bool Canvas::indicesInRange(double x1, double x2, size_t& i1, size_t& i2)
{
    // ensure x1 <= x2
    if (x1 > x2) std::swap(x1, x2);
    // get the raw indices based on pixel center; these may be out of range (including negative)
    long index1 = static_cast<long>(_n * x1 + 0.5);
    long index2 = static_cast<long>(_n * x2 - 0.5);
    // if the range is empty (because narrower than the pixel size), use the intersected pixel
    if (index2 < index1) index1 = index2 = static_cast<long>(_n * 0.5 * (x1 + x2));
    // if the range is outside of the canvas, report failure
    if (index2 < 0 || index1 >= static_cast<long>(_n)) return false;
    // otherwise, clip the indices to the canvas
    if (index1 < 0) index1 = 0;
    if (index2 >= static_cast<long>(_n)) index2 = static_cast<long>(_n) - 1;
    i1 = static_cast<size_t>(index1);
    i2 = static_cast<size_t>(index2);
    return true;
}

////////////////////////////////////////////////////////////////////

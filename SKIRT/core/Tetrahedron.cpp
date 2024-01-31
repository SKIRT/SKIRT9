#include "Tetrahedron.hpp"
#include "FatalError.hpp"
#include <iostream>

Tetra::Tetra(const std::array<Vec*, 4>& vertices, const std::array<Face, 4>& faces) : _vertices(vertices), _faces(faces)
{
    double xmin = DBL_MAX;
    double ymin = DBL_MAX;
    double zmin = DBL_MAX;
    double xmax = -DBL_MAX;
    double ymax = -DBL_MAX;
    double zmax = -DBL_MAX;
    for (const Vec* vertex : _vertices)
    {
        xmin = min(xmin, vertex->x());
        ymin = min(ymin, vertex->y());
        zmin = min(zmin, vertex->z());
        xmax = max(xmax, vertex->x());
        ymax = max(ymax, vertex->y());
        zmax = max(zmax, vertex->z());
    }
    setExtent(Box(xmin, ymin, zmin, xmax, ymax, zmax));

    _volume = 1 / 6.
              * abs(Vec::dot(Vec::cross(*_vertices[1] - *_vertices[0], *_vertices[2] - *_vertices[0]),
                             *_vertices[3] - *_vertices[0]));

    for (int i = 0; i < 4; i++) _centroid += *_vertices[i];
    _centroid /= 4;

    for (int f = 0; f < 4; f++)
    {
        std::array<int, 3> cv = clockwiseVertices(f);
        Vec e12 = *vertices[cv[1]] - *vertices[cv[0]];
        Vec e13 = *vertices[cv[2]] - *vertices[cv[0]];
        Vec normal = Vec::cross(e12, e13);
        normal /= normal.norm();

        Face& face = _faces[f];
        face._normal = normal;
    }

    // this convention makes edges go clockwise around leaving rays from inside the tetrahedron
    // so their plucker products are all positive if the ray leaves
    const Vec e01 = *_vertices[1] - *_vertices[0];
    const Vec e02 = *_vertices[2] - *_vertices[0];
    const Vec e03 = *_vertices[3] - *_vertices[0];
    double orientation = Vec::dot(Vec::cross(e01, e02), e03);
    if (orientation < 0)
    {
        printf("ORIENTATION SWITCHED!!!!!!!!!");
        // swap last 2, this means first 2 indices can be ordered i < j
    }
}

////////////////////////////////////////////////////////////////////

Vec Tetra::getEdge(int t1, int t2) const
{
    return *_vertices[t2] - *_vertices[t1];
}

////////////////////////////////////////////////////////////////////

bool Tetra::inside(const Position& bfr) const
{
    if (!Box::contains(bfr)) return false;

    /*
    face: normals for which the other vertex has a positive dot product with
    3:*02 x 01*| 10 x 12 | 21 x 20
    2: 13 x 10 |*01 x 03*| 30 x 31
    1: 20 x 23 | 32 x 30 |*03 x 02*
    0: 31 x 32 | 23 x 21 |*12 x 13* // last one doesn't matter
    */

    // optimized version (probably not that much better)
    Vec e0p = bfr - *_vertices[0];
    Vec e02 = getEdge(0, 2);
    Vec e01 = getEdge(0, 1);
    if (Vec::dot(Vec::cross(e02, e01), e0p) > 0)  // 02 x 01
        return false;

    Vec e03 = *_vertices[3] - *_vertices[0];
    if (Vec::dot(Vec::cross(e01, e03), e0p) > 0)  // 01 x 03
        return false;

    if (Vec::dot(Vec::cross(e03, e02), e0p) > 0)  // 03 x 02
        return false;

    Vec e1p = bfr - *_vertices[1];
    Vec e12 = getEdge(1, 2);
    Vec e13 = getEdge(1, 3);
    return Vec::dot(Vec::cross(e12, e13), e1p) < 0;  // 12 x 13

    // checks 3 edges too many but very simple
    // for (int face = 0; face < 4; face++)
    // {
    //     std::array<int, 3> t = clockwiseVertices(face);
    //     Vec& v0 = *_vertices[t[0]];
    //     Vec& clock = *_vertices[t[1]];
    //     Vec& counter = *_vertices[t[2]];
    //     Vec normal = Vec::cross(counter - v0, clock - v0);
    //     if (Vec::dot(normal, bfr - v0) < 0)  // is pos on the same side as v3
    //     {
    //         return false;
    //     }
    // }
    // return true;
}

////////////////////////////////////////////////////////////////////

double Tetra::volume() const
{
    return _volume;
}

////////////////////////////////////////////////////////////////////

const Array& Tetra::properties()
{
    return _properties;
}

////////////////////////////////////////////////////////////////////

double Tetra::generateBarycentric(double& s, double& t, double& u)
{
    // https://vcg.isti.cnr.it/activities/OLD/geometryegraphics/pointintetraedro.html
    if (s + t > 1.0)
    {  // cut'n fold the cube into a prism

        s = 1.0 - s;
        t = 1.0 - t;
    }
    if (t + u > 1.0)
    {  // cut'n fold the prism into a tetrahedron

        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
    }
    else if (s + t + u > 1.0)
    {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
    }
    return 1 - u - t - s;
}

////////////////////////////////////////////////////////////////////

Position Tetra::generatePosition(double s, double t, double u) const
{
    double w = Tetra::generateBarycentric(s, t, u);

    return Position(w * *_vertices[0] + u * *_vertices[1] + t * *_vertices[2] + s * *_vertices[3]);
}

////////////////////////////////////////////////////////////////////

std::array<int, 3> Tetra::clockwiseVertices(int face)
{
    std::array<int, 3> cv = {(face + 3) % 4, (face + 2) % 4, (face + 1) % 4};
    // if face is even we should swap two edges
    if (face % 2 == 0) std::swap(cv[0], cv[2]);
    return cv;
}

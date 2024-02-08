#include "Tetrahedron.hpp"
#include "FatalError.hpp"
#include <iostream>

Tetra::Tetra(const std::array<Vec*, 4>& vertices, const std::array<Face, 4>& faces, const Array& prop)
    : Tetra(vertices, faces)
{
    _properties = prop;
}

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

    // volume
    _volume = 1 / 6.
              * abs(Vec::dot(Vec::cross(*_vertices[1] - *_vertices[0], *_vertices[2] - *_vertices[0]),
                             *_vertices[3] - *_vertices[0]));

    // barycenter
    for (int i = 0; i < 4; i++) _centroid += *_vertices[i];
    _centroid /= 4;

    // calculate normal facing out
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
}

////////////////////////////////////////////////////////////////////

Vec Tetra::getEdge(int t1, int t2) const
{
    return *_vertices[t2] - *_vertices[t1];
}

////////////////////////////////////////////////////////////////////

bool Tetra::inside(const Position& bfr) const
{
    // since we use the k-d tree this will only slow the CellIndex
    // if (!Box::contains(bfr)) return false;

    // could optimize this slightly by using same vertex for 3 faces and do final face seperately

    for (int f = 0; f < 4; f++)
    {
        const Face& face = _faces[f];
        const Vec* vertex = _vertices[(f + 1) % 4];  // any vertex that is on the face

        // if point->face is opposite direction as the outward pointing normal, the point is outside
        if (Vec::dot(*vertex - bfr, face._normal) < 0) return false;
    }
    return true;
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

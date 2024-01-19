#include "Tetrahedron.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

Plucker::Plucker() {}

////////////////////////////////////////////////////////////////////

Plucker::Plucker(const Vec& pos, const Vec& dir) : U(dir), V(Vec::cross(dir, pos)) {}

////////////////////////////////////////////////////////////////////

inline double Plucker::dot(const Plucker& a, const Plucker& b)
{
    return Vec::dot(a.U, b.V) + Vec::dot(b.U, a.V);
}

////////////////////////////////////////////////////////////////////

Edge::Edge(int i1, int i2, const Vec* v1, const Vec* v2) : Plucker(*v1, *v2 - *v1), i1(i1), i2(i2) {}

////////////////////////////////////////////////////////////////////

Tetra::Tetra(const std::array<Vec*, 4>& vertices, const std::array<int, 4>& indices,
             const std::array<int, 4>& neighbors, const std::array<Edge*, 6>& edges)
    : _vertices(vertices), _indices(indices), _edges(edges), _neighbors(neighbors)
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
}

////////////////////////////////////////////////////////////////////

double Tetra::getProd(const Plucker& ray, int t1, int t2) const
{
    int e = (std::min(t1, t2) == 0) ? std::max(t1, t2) - 1 : t1 + t2;
    Edge* edge = _edges[e];

    // not the same order -> *-1
    // beginning of t1 == beginning of edge (= same order)
    return (_indices[t1] == edge->i1 ? 1 : -1) * Plucker::dot(ray, *edge);
}

////////////////////////////////////////////////////////////////////

bool Tetra::intersects(std::array<double, 3>& barycoords, const Plucker& ray, int face, bool leaving) const
{
    std::array<int, 3> t = clockwiseVertices(face);

    double sum = 0;
    for (int i = 0; i < 3; i++)
    {
        // edges: 12, 20, 01
        // verts:  0,  1,  2
        double prod = getProd(ray, t[(i + 1) % 3], t[(i + 2) % 3]);
        if (leaving != (prod <= 0)) return false;  // change this so for both leavig and entering prod=0 works
        barycoords[i] = prod;
        sum += prod;
    }
    if (sum != 0)
        for (int i = 0; i < 3; i++) barycoords[i] /= sum;

    return true;
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
    Vec e02 = *_vertices[2] - *_vertices[0];
    Vec e01 = *_vertices[1] - *_vertices[0];
    if (Vec::dot(Vec::cross(e02, e01), e0p) > 0)  // 02 x 01
        return false;

    Vec e03 = *_vertices[3] - *_vertices[0];
    if (Vec::dot(Vec::cross(e01, e03), e0p) > 0)  // 01 x 03
        return false;

    if (Vec::dot(Vec::cross(e03, e02), e0p) > 0)  // 03 x 02
        return false;

    Vec e1p = bfr - *_vertices[1];
    Vec e12 = *_vertices[2] - *_vertices[1];
    Vec e13 = *_vertices[3] - *_vertices[1];
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

Vec Tetra::calcExit(const std::array<double, 3>& barycoords, int face) const
{
    std::array<int, 3> t = Tetra::clockwiseVertices(face);
    Vec exit;
    for (int i = 0; i < 3; i++) exit += *_vertices[t[i]] * barycoords[i];
    return exit;
}

////////////////////////////////////////////////////////////////////

Vec Tetra::centroid() const
{
    Vec pos;
    for (int i = 0; i < 4; i++) pos += *_vertices[i];
    pos /= 4;
    return pos;
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
    std::array<int, 3> t = {(face + 3) % 4, (face + 2) % 4, (face + 1) % 4};
    // if face is even we should swap two edges
    if (face % 2 == 0) std::swap(t[0], t[2]);
    return t;
}

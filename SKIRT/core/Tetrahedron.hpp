#include "Array.hpp"
#include "Box.hpp"
#include "Position.hpp"

class Plucker
{
private:
    Vec U, V;

public:
    Plucker();

    Plucker(const Vec& pos, const Vec& dir);

    // permuted inner product
    static inline double dot(const Plucker& a, const Plucker& b);
};

class Edge : public Plucker
{
public:
    const int i1, i2;
    Edge(int i1, int i2, const Vec* v1, const Vec* v2);
};

class Tetra : public Box
{
private:
    double _volume;
    Array _properties;

public:
    std::array<Vec*, 4> _vertices;
    std::array<int, 4> _indices;
    std::array<Edge*, 6> _edges;
    std::array<int, 4> _neighbors;

public:
    Tetra(const std::array<Vec*, 4>& vertices, const std::array<int, 4>& indices, const std::array<int, 4>& neighbors,
          const std::array<Edge*, 6>& edges);

    Tetra(Vec* va, Vec* vb, Vec* vc, Vec* vd);

    double getProd(const Plucker& ray, int t1, int t2) const;

    bool intersects(std::array<double, 3>& prods, const Plucker& ray, int face, bool leaving = true) const;

    bool inside(const Position& bfr) const;

    Vec calcExit(const std::array<double, 3>& barycoords, int face) const;

    Vec centroid() const;

    double volume() const;

    const Array& properties();

    Position generatePosition(double s, double t, double u) const;

    static std::array<int, 3> clockwiseVertices(int face);
};

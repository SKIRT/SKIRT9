#include "Position.hpp"
#include "Array.hpp"
#include "Box.hpp"
#include <array>

struct Face
{
    Face() {}

    // normals are calculated in the constructor of Tetra
    Face(int ntetra, int nface) : _ntetra(ntetra), _nface(nface) {}

    Vec _normal;
    int _ntetra;
    int _nface;
};

class Tetra : public Box
{
private:
    double _volume;
    Array _properties;

public:
    const std::array<Vec*, 4> _vertices;
    std::array<Face, 4> _faces;
    Vec _centroid;

public:
    Tetra(const std::array<Vec*, 4>& vertices, const std::array<Face, 4>& faces, const Array& prop);

    Tetra(const std::array<Vec*, 4>& vertices, const std::array<Face, 4>& faces);

    bool inside(const Position& bfr) const;

    Vec getEdge(int t1, int t2) const;

    double volume() const;

    const Array& properties();

    static double generateBarycentric(double& s, double& t, double& u);

    Position generatePosition(double s, double t, double u) const;

    /**
     * @brief This gives the clockwise vertices of a given face looking from inside the tetrahedron.
     * The Plücker products are thus all negative for a leaving ray and positive for an entering ray from outside
     * 
     * @param face 
     * @return std::array<int, 3> 
     */
    static std::array<int, 3> clockwiseVertices(int face);
};

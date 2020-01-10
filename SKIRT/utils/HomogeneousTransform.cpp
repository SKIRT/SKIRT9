/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HomogeneousTransform.hpp"

//////////////////////////////////////////////////////////////////////

HomogeneousTransform::HomogeneousTransform()
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) M[i][j] = i == j ? 1. : 0.;
}

void HomogeneousTransform::translate(double x, double y, double z)
{
    HomogeneousTransform transform;
    transform.M[3][0] = x;
    transform.M[3][1] = y;
    transform.M[3][2] = z;
    concatenate(transform);
}

//////////////////////////////////////////////////////////////////////

void HomogeneousTransform::scale(double x, double y, double z)
{
    HomogeneousTransform transform;
    transform.M[0][0] = x;
    transform.M[1][1] = y;
    transform.M[2][2] = z;
    concatenate(transform);
}

//////////////////////////////////////////////////////////////////////

void HomogeneousTransform::rotateX(double cos, double sin)
{
    HomogeneousTransform transform;
    transform.M[1][1] = cos;
    transform.M[2][2] = cos;
    transform.M[1][2] = -sin;
    transform.M[2][1] = sin;
    concatenate(transform);
}

//////////////////////////////////////////////////////////////////////

void HomogeneousTransform::rotateY(double cos, double sin)
{
    HomogeneousTransform transform;
    transform.M[0][0] = cos;
    transform.M[2][2] = cos;
    transform.M[0][2] = -sin;
    transform.M[2][0] = sin;
    concatenate(transform);
}

//////////////////////////////////////////////////////////////////////

void HomogeneousTransform::rotateZ(double cos, double sin)
{
    HomogeneousTransform transform;
    transform.M[0][0] = cos;
    transform.M[1][1] = cos;
    transform.M[0][1] = -sin;
    transform.M[1][0] = sin;
    concatenate(transform);
}

//////////////////////////////////////////////////////////////////////

void HomogeneousTransform::perspectiveZ(double f)
{
    HomogeneousTransform transform;
    transform.M[2][3] = 1. / f;
    transform.M[3][2] = -f;
    transform.M[3][3] = 0.;
    concatenate(transform);
}

//////////////////////////////////////////////////////////////////////

void HomogeneousTransform::concatenate(const HomogeneousTransform& transform)
{
    HomogeneousTransform copyofthis(*this);
    typedef double matrix[4][4];
    const matrix& M1 = copyofthis.M;
    const matrix& M2 = transform.M;

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            M[i][j] = M1[i][0] * M2[0][j] + M1[i][1] * M2[1][j] + M1[i][2] * M2[2][j] + M1[i][3] * M2[3][j];
}

//////////////////////////////////////////////////////////////////////

void HomogeneousTransform::transform(double x, double y, double z, double w, double& outx, double& outy, double& outz,
                                     double& outw) const
{
    outx = x * M[0][0] + y * M[1][0] + z * M[2][0] + w * M[3][0];
    outy = x * M[0][1] + y * M[1][1] + z * M[2][1] + w * M[3][1];
    outz = x * M[0][2] + y * M[1][2] + z * M[2][2] + w * M[3][2];
    outw = x * M[0][3] + y * M[1][3] + z * M[2][3] + w * M[3][3];
}

//////////////////////////////////////////////////////////////////////

Vec HomogeneousTransform::transform(Vec p) const
{
    double x, y, z, w;
    transform(p.x(), p.y(), p.z(), 1., x, y, z, w);
    return Vec(x / w, y / w, z / w);
}

//////////////////////////////////////////////////////////////////////

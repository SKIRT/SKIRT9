/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef QUADRICS_HPP
#define QUADRICS_HPP

#include "Vec.hpp"

////////////////////////////////////////////////////////////////////

/** This static class offers a collection of low-level geometric building blocks: raising a value
    to the second or third power, solving a quadratic equation, and a set of elementary
    sphere/plane tests and ray/quadric intersections centered on or relative to the origin.

    All implementations are provided inline in the header. */
class Quadrics final
{
    // small value
    constexpr static double EPS = 1e-12;

    //======================== Powers =======================

public:
    /** This function returns \f$x^2\f$. */
    static double square(double x) { return x * x; }

    /** This function returns \f$x^3\f$. */
    static double pow3(double x) { return x * x * x; }

    /** This function returns \f$x_1^3 - x_0^3 = (x_1-x_0)(x_1^2 + x_1 x_0 + x_0^2)\f$. The second
        form is used because it is more numerically stable. */
    static double pow3(double x0, double x1) { return (x1 - x0) * (x1 * x1 + x1 * x0 + x0 * x0); }

    //======================== Quadratic equations =======================

public:
    /** This function determines the solutions of \f$x^2 + 2bx + c = 0\f$. If there are two
        distinct real solutions, they are stored in the arguments x1 and x2, and the function
        returns 2. Otherwise, i.e. if there are no real solutions or there is just one real
        solution, x1 and x2 remain unchanged, and the function returns 0. */
    static int distinctSolutions(double b, double c, double& x1, double& x2)
    {
        if (b * b > c)  // if discriminant is strictly positive, there are two distinct real solutions
        {
            if (b > 0)  // x1 is always negative
            {
                x1 = -b - sqrt(b * b - c);
                x2 = c / x1;
            }
            else  // x2 is always positive
            {
                x2 = -b + sqrt(b * b - c);
                x1 = c / x2;
            }
            return 2;
        }
        return 0;
    }

    /** This function determines the solutions of \f$ax^2 + 2bx + c = 0\f$. If there are two
        distinct real solutions, they are stored in the arguments x1 and x2, and the function
        returns 2. If the equation degenerates to a linear equation with a single finite solution,
        this solution is stored in the argument x1. In that case, x2 remains unchanged and the
        function returns 1. Otherwise, i.e. if there are no real solutions or the non-generate
        quadratic equation has two identical solutions, x1 and x2 remain unchanged, and the
        function returns 0. */
    static int distinctSolutions(double a, double b, double c, double& x1, double& x2)
    {
        if (std::abs(a) > EPS) return distinctSolutions(b / a, c / a, x1, x2);
        double x = -0.5 * c / b;
        if (std::isfinite(x))
        {
            x1 = x;
            return 1;
        }
        return 0;
    }

    /** This function returns the smallest positive solution of \f$x^2 + 2bx + c = 0\f$, or zero if
        there is no positive solution. */
    static double smallestPositiveSolution(double b, double c)
    {
        if (b * b > c)  // if discriminant is negative, there are no real solutions
        {
            if (b > 0.)  // x1 is always negative; x2 is positive only if c<0
            {
                if (c < 0.)
                {
                    double x1 = -b - sqrt(b * b - c);
                    return c / x1;
                }
            }
            else  // x2 is always positive; x1 is positive only if c>0
            {
                double x2 = -b + sqrt(b * b - c);
                if (c > 0.)
                {
                    double x1 = c / x2;
                    if (x1 < x2) return x1;
                }
                return x2;
            }
        }
        return 0.;
    }

    /** This function returns the smallest positive solution of \f$ax^2 + 2bx + c = 0\f$, or zero
        if there is no positive solution. If the equation degenerates to a linear equation, this
        equation is solved instead. */
    static double smallestPositiveSolution(double a, double b, double c)
    {
        if (std::abs(a) > EPS) return smallestPositiveSolution(b / a, c / a);
        double x = -0.5 * c / b;
        if (std::isfinite(x) && x > 0.) return x;
        return 0.;
    }

    //======================== Spheres =======================

public:
    /** This function returns the volume of a sphere with given radius. */
    static double volumeSphere(double radius) { return (4. / 3.) * M_PI * pow3(radius); }

    /** This function returns true if the given position is inside the sphere with given center and
        radius, and false otherwise (a position exactly on the sphere's surface counts as inside).
        */
    static bool isPositionInSphere(Vec p, Vec center, double radius) { return (p - center).norm2() <= square(radius); }

    /** This function returns true if the specified spheres overlap or touch, and false otherwise.
        */
    static bool doSpheresOverlap(Vec center1, double radius1, Vec center2, double radius2)
    {
        return (center1 - center2).norm2() <= square(radius1 + radius2);
    }

    /** This function returns true if the sphere with given center and radius is fully inside the
        spherical shell, centered on the origin, with given inner and outer radii (no straggling
        nor touching allowed). */
    static bool isSphereInShell(Vec center, double radius, double rmin, double rmax)
    {
        const double rho = center.norm();
        return rho - radius > rmin && rho + radius < rmax;
    }

    /** This function determines the circle formed by intersecting a sphere (given its radius and
        its center's coordinate along the axis perpendicular to the plane) with a coordinate plane
        through the origin. It returns true and sets circleRadius if the sphere reaches the plane;
        it returns false (leaving circleRadius unchanged) if it doesn't. The circle's center in the
        plane is simply the sphere center's own two in-plane coordinates, so the caller already has
        those and doesn't need them returned here. */
    static bool sphereIntersectsPlane(double sphereRadius, double outOfPlaneCenterCoord, double& circleRadius)
    {
        double r2 = square(sphereRadius) - square(outOfPlaneCenterCoord);
        if (r2 <= 0.) return false;
        circleRadius = std::sqrt(r2);
        return true;
    }

    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        {\bf{k}})\f$ -- with \f${\bf{k}}\f$ assumed to be a unit vector -- with the sphere of given
        radius centered at the origin, or zero if there is no intersection ahead along the ray
        (either no real intersection, or both intersections lie behind the ray's origin). */
    static double firstIntersectionSphere(Vec r0, Vec k, double radius)
    {
        return smallestPositiveSolution(Vec::dot(r0, k), r0.norm2() - square(radius));
    }

    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        {\bf{k}})\f$ with the sphere of given center and radius, or zero if there is no
        intersection ahead along the ray. This is the arbitrary-center counterpart of the other
        firstIntersectionSphere() overload, obtained by translating the ray so that the sphere's
        center coincides with the origin. */
    static double firstIntersectionSphere(Vec r0, Vec k, Vec center, double radius)
    {
        return firstIntersectionSphere(r0 - center, k, radius);
    }

    //======================== Cones and meridional planes through the origin =======================

public:
    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        {\bf{k}})\f$ with the double cone, centered on the origin and aligned with the z-axis, that
        has the given cosine of its opening angle, or zero if there is no intersection ahead along
        the ray. The degenerate cone with zero cosine (i.e. the xy-plane) is handled as a special
        case. */
    static double firstIntersectionCone(Vec r0, Vec k, double cosTheta)
    {
        double c = cosTheta;
        return c ? smallestPositiveSolution(c * c - k.z() * k.z(), c * c * Vec::dot(r0, k) - r0.z() * k.z(),
                                            c * c * r0.norm2() - r0.z() * r0.z())
                 : -r0.z() / k.z();  // degenerate cone identical to xy-plane
    }

    /** This function returns the distance to the intersection of the ray \f$({\bf{r}}_0,
        {\bf{k}})\f$ with the meridional half-plane, through the origin and the z-axis, at the
        azimuth with given sine and cosine, or zero if the ray direction is parallel to the plane.
        */
    static double intersectionMeridionalPlane(Vec r0, Vec k, double sinPhi, double cosPhi)
    {
        double q = k.x() * sinPhi - k.y() * cosPhi;
        if (std::abs(q) < EPS) return 0.;
        return -(r0.x() * sinPhi - r0.y() * cosPhi) / q;
    }

    //======================== Cylinders around the z-axis and horizontal planes =======================

public:
    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        {\bf{k}})\f$ with the infinite cylinder of given radius, centered on and aligned with the
        z-axis, or zero if there is no intersection ahead along the ray (either no real
        intersection, both intersections lie behind the ray's origin, or the ray direction is
        parallel to the z-axis). The squared length \f$k_x^2+k_y^2\f$ of the projection of
        \f${\bf{k}}\f$ onto the xy-plane is passed in as kq2 rather than recomputed, because the
        caller typically needs the first intersection with several cylinders of different radii for
        the same ray, all sharing the same kq2. */
    static double firstIntersectionCylinder(Vec r0, Vec k, double kq2, double radius)
    {
        if (std::abs(kq2) < EPS) return 0.;
        double b = r0.x() * k.x() + r0.y() * k.y();
        double c = square(r0.x()) + square(r0.y()) - square(radius);
        return smallestPositiveSolution(b / kq2, c / kq2);
    }

    /** This function returns the distance to the intersection of the ray \f$({\bf{r}}_0,
        {\bf{k}})\f$ with the horizontal plane, perpendicular to the z-axis, at the given z
        coordinate, or zero if the ray direction is parallel to the plane. */
    static double intersectionHorizontalPlane(Vec r0, Vec k, double zPlane)
    {
        if (std::abs(k.z()) < EPS) return 0.;
        return (zPlane - r0.z()) / k.z();
    }
};

////////////////////////////////////////////////////////////////////

#endif

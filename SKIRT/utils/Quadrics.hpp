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
    static double cube(double x) { return x * x * x; }

    /** This function returns \f$x_1^3 - x_0^3 = (x_1-x_0)(x_1^2 + x_1 x_0 + x_0^2)\f$. The second
        form is used because it is more numerically stable. */
    static double cube(double x0, double x1) { return (x1 - x0) * (x1 * x1 + x1 * x0 + x0 * x0); }

    /** This function returns \f$\sqrt{a^2-b^2}\f$ for \f$|a|\ge|b|\f$, computed as
        \f$\sqrt{(a-b)(a+b)}\f$ rather than directly from \f$a^2-b^2\f$, to avoid the loss of
        precision that direct evaluation would incur through catastrophic cancellation when
        \f$|a|\f$ and \f$|b|\f$ are close. The argument to the square root is clamped to zero in
        case rounding error would otherwise make it (slightly) negative even though the exact
        value is nonnegative. */
    static double sqrtDiffSquares(double a, double b) { return sqrt(max(0., (a - b) * (a + b))); }

    //======================== Quadratic equations =======================

public:
    /** This function determines the solutions of \f$x^2 + 2bx + c = 0\f$. If there are two
        distinct real solutions, they are stored in the arguments x1 and x2, and the function
        returns 2. Otherwise, i.e. if there are no real solutions or there is just one real
        solution, x1 and x2 remain unchanged, and the function returns 0.

        Evaluating \f$-b \pm \sqrt{b^2-c}\f$ directly for both roots would lose precision through
        catastrophic cancellation whenever the two terms nearly cancel for one of the roots. To
        avoid this, only the root for which \f$-b\f$ and the square root have the same sign --
        so that they add rather than cancel -- is evaluated directly: \f$-b-\sqrt{b^2-c}\f$ if
        \f$b>0\f$ (since \f$-b\f$ is then negative), or \f$-b+\sqrt{b^2-c}\f$ otherwise. The other
        root is then obtained from Vieta's formula \f$x_1 x_2 = c\f$. */
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
        there is no positive solution.

        As in distinctSolutions(), only the root for which \f$-b\f$ and \f$\sqrt{b^2-c}\f$ have
        the same sign is evaluated directly -- \f$-b-\sqrt{b^2-c}\f$ if \f$b>0\f$, or
        \f$-b+\sqrt{b^2-c}\f$ if \f$b<0\f$ -- avoiding the catastrophic cancellation that
        evaluating both roots directly would incur; the other root, when it is needed, is
        obtained from it through Vieta's formula \f$x_1 x_2 = c\f$ rather than by evaluating the
        numerically unsafe combination directly. */
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
    /** This function returns the volume of a sphere with given radius, \f$V = \frac{4}{3}\pi
        r^3\f$. */
    static double volumeSphere(double radius) { return (4. / 3.) * M_PI * cube(radius); }

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

private:
    // returns the minimum, over t ranging across the interval [tlo,thi], of amp*cos(t-phase),
    // where amp is assumed nonnegative; used by isSphericalCellInSphere() to locate the point on a bounded
    // arc (in azimuth or in polar angle) that lies farthest, in the sense of this cosine
    // projection, from a given reference direction
    static double minCosineOverInterval(double amp, double phase, double tlo, double thi)
    {
        double antipodal = phase + M_PI;                        // phase is in (-pi,pi], so this lies in (0,2pi]
        if (antipodal > M_PI) antipodal -= 2. * M_PI;           // wrap back into (-pi,pi]
        if (antipodal >= tlo && antipodal <= thi) return -amp;  // the global minimum falls inside the arc
        return min(amp * cos(tlo - phase), amp * cos(thi - phase));
    }

public:
    /** This function returns true if the "spherical box" cell -- the region with radius between
        rmin and rmax, polar angle (colatitude) between thetaMin and thetaMax, and azimuth between
        phiMin and phiMax, all measured around the origin -- lies fully inside the sphere with the
        given center and radius (touching allowed).

        The radial extent is handled exactly: for a fixed direction, the squared distance to a
        point at radius r is a convex function of r, so its maximum over [rmin,rmax] is always
        attained at one of the two endpoints, regardless of curvature. What remains is finding the
        direction within the cell's angular box that lies farthest, in projection, from the
        sphere's center. Because the azimuthal term of that projection is scaled by sin(theta),
        which is nonnegative throughout the box, the farthest azimuth can be found first,
        independently of theta, and the result then feeds a second, one-dimensional search over
        theta. Each of those two searches reduces to checking whether the antipodal direction falls
        within the relevant bound, and otherwise evaluating the two boundary values -- never a full
        numerical optimization. */
    static bool isSphericalCellInSphere(double rmin, double rmax, double thetaMin, double thetaMax, double phiMin,
                                        double phiMax, Vec center, double radius)
    {
        double phiC = atan2(center.y(), center.x());
        double Cxy = sqrt(square(center.x()) + square(center.y()));
        double farthestAzimuthProjection = minCosineOverInterval(Cxy, phiC, phiMin, phiMax);

        double psi = atan2(farthestAzimuthProjection, center.z());
        double R = sqrt(square(farthestAzimuthProjection) + square(center.z()));
        double farthestProjection = minCosineOverInterval(R, psi, thetaMin, thetaMax);

        double c2 = center.norm2();
        double d2min = square(rmin) - 2. * rmin * farthestProjection + c2;
        double d2max = square(rmax) - 2. * rmax * farthestProjection + c2;
        return max(d2min, d2max) <= square(radius);
    }

    /** This function determines the circle formed by intersecting a sphere (given its radius and
        its center's coordinate along the axis perpendicular to the plane) with a coordinate plane
        through the origin. It returns true and sets circleRadius if the sphere reaches the plane;
        it returns false (leaving circleRadius unchanged) if it doesn't. The circle's center in the
        plane is simply the sphere center's own two in-plane coordinates, so the caller already has
        those and doesn't need them returned here.

        By the Pythagorean relation between the sphere's radius \f$R\f$, its center's
        out-of-plane offset \f$d\f$, and the radius \f$\rho\f$ of the circle cut out of the plane,
        \f[ \rho = \sqrt{R^2 - d^2}, \f] which is real, and thus returned, only if \f$|d|<R\f$. */
    static bool sphereIntersectsPlane(double sphereRadius, double outOfPlaneCenterCoord, double& circleRadius)
    {
        if (std::abs(outOfPlaneCenterCoord) >= sphereRadius) return false;
        circleRadius = sqrtDiffSquares(sphereRadius, outOfPlaneCenterCoord);
        return true;
    }

    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the sphere of given radius centered at the origin, or zero if there is no
        intersection ahead along the ray (either no real intersection, or both intersections lie
        behind the ray's origin).

        Substituting \f${\bf{p}}(s) = {\bf{r}}_0 + s\hat{\bf{k}}\f$ into the sphere's defining
        equation \f$|{\bf{p}}|^2 = \text{radius}^2\f$ gives the quadratic \f[ s^2 + 2s\,
        ({\bf{r}}_0\cdot\hat{\bf{k}}) + (r_0^2 - \text{radius}^2) = 0, \f] whose smallest positive
        solution, if any, is the returned distance. */
    static double firstIntersectionSphere(Vec r0, Vec k, double radius)
    {
        return smallestPositiveSolution(Vec::dot(r0, k), r0.norm2() - square(radius));
    }

    /** This function determines the solutions of the intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the sphere described in firstIntersectionSphere(), i.e. the two-root
        counterpart of that function: it solves the same defining equation and quadratic, but
        through distinctSolutions() rather than smallestPositiveSolution(), so that both roots are
        returned rather than only the smallest positive one. The return value has the same meaning
        as for distinctSolutions(): 2 for two distinct solutions (stored in x1 and x2), or 0 if
        there are none. */
    static int distinctIntersectionsSphere(Vec r0, Vec k, double radius, double& x1, double& x2)
    {
        return distinctSolutions(Vec::dot(r0, k), r0.norm2() - square(radius), x1, x2);
    }

    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the sphere of given center and radius, or zero if there is no
        intersection ahead along the ray. This is the arbitrary-center counterpart of the other
        firstIntersectionSphere() overload, obtained by translating the ray so that the sphere's
        center coincides with the origin -- i.e. by substituting \f${\bf{r}}_0 - \text{center}\f$
        for \f${\bf{r}}_0\f$ in that overload's defining equation. */
    static double firstIntersectionSphere(Vec r0, Vec k, Vec center, double radius)
    {
        return firstIntersectionSphere(r0 - center, k, radius);
    }

    //======================== Cones and meridional planes through the origin =======================

public:
    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the double cone, centered on the origin and aligned with the z-axis, that
        has the given cosine of its opening angle, or zero if there is no intersection ahead along
        the ray. The degenerate cone with zero cosine (i.e. the xy-plane) is handled as a special
        case.

        Substituting \f${\bf{p}}(s) = {\bf{r}}_0 + s\hat{\bf{k}}\f$ into the cone's defining
        equation \f$\cos^2\theta\,|{\bf{p}}|^2 = p_z^2\f$ gives the quadratic \f[ (\cos^2\theta -
        k_z^2)\,s^2 + 2s\,(\cos^2\theta\,({\bf{r}}_0\cdot\hat{\bf{k}}) - r_{0,z}\,k_z) +
        (\cos^2\theta\,r_0^2 - r_{0,z}^2) = 0, \f] whose smallest positive solution, if any, is the
        returned distance. */
    static double firstIntersectionCone(Vec r0, Vec k, double cosTheta)
    {
        double c = cosTheta;
        if (c)
            return smallestPositiveSolution(c * c - k.z() * k.z(), c * c * Vec::dot(r0, k) - r0.z() * k.z(),
                                            c * c * r0.norm2() - r0.z() * r0.z());

        // degenerate cone identical to xy-plane; intersection exists only if not parallel to it
        if (std::abs(k.z()) < EPS) return 0.;
        double s = -r0.z() / k.z();
        return s > 0. ? s : 0.;
    }

    /** This function determines the solutions of the intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the double cone described in firstIntersectionCone(), i.e. the
        two-root counterpart of that function: it solves the same defining equation and quadratic,
        but through distinctSolutions() rather than smallestPositiveSolution(), so that both roots
        (or the single root of a linearly-degenerate case) are returned rather than only the
        smallest positive one. The degenerate cone with zero cosine (i.e. the xy-plane) is handled
        as a special case, with at most one solution, obtained only if the ray is not parallel to
        it. The return value has the same meaning as for distinctSolutions(): 2 for two distinct
        solutions (stored in x1 and x2), 1 for a single solution (stored in x1), or 0 for none. */
    static int distinctIntersectionsCone(Vec r0, Vec k, double cosTheta, double& x1, double& x2)
    {
        double c2 = square(cosTheta);
        if (!c2)
        {
            // degenerate cone identical to xy-plane; at most one solution, and only if not parallel
            if (std::abs(k.z()) < EPS) return 0;
            x1 = -r0.z() / k.z();
            return 1;
        }
        return distinctSolutions(c2 - square(k.z()), c2 * Vec::dot(r0, k) - r0.z() * k.z(),
                                 c2 * r0.norm2() - square(r0.z()), x1, x2);
    }

    /** This function returns the distance to the intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the meridional half-plane, through the origin and the z-axis, at the
        azimuth with given sine and cosine, or zero if the ray direction is parallel to the plane.

        Substituting \f${\bf{p}}(s) = {\bf{r}}_0 + s\hat{\bf{k}}\f$ into the plane's defining
        equation \f$p_x\sin\varphi - p_y\cos\varphi = 0\f$ and solving for \f$s\f$ gives \f[ s =
        -\frac{r_{0,x}\sin\varphi - r_{0,y}\cos\varphi}{k_x\sin\varphi - k_y\cos\varphi}. \f] */
    static double intersectionMeridionalPlane(Vec r0, Vec k, double sinPhi, double cosPhi)
    {
        double q = k.x() * sinPhi - k.y() * cosPhi;
        if (std::abs(q) < EPS) return 0.;
        return -(r0.x() * sinPhi - r0.y() * cosPhi) / q;
    }

    //======================== Cylinders around the z-axis and horizontal planes =======================

public:
    /** This function returns the distance to the first intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the infinite cylinder of given radius, centered on and aligned with the
        z-axis, or zero if there is no intersection ahead along the ray (either no real
        intersection, both intersections lie behind the ray's origin, or the ray direction is
        parallel to the z-axis). The squared length \f$k_x^2+k_y^2\f$ of the projection of
        \f$\hat{\bf{k}}\f$ onto the xy-plane is passed in as kq2 rather than recomputed, because the
        caller typically needs the first intersection with several cylinders of different radii for
        the same ray, all sharing the same kq2.

        Substituting \f${\bf{p}}(s) = {\bf{r}}_0 + s\hat{\bf{k}}\f$ into the cylinder's defining
        equation \f$p_x^2+p_y^2 = \text{radius}^2\f$ gives the quadratic \f[ (k_x^2+k_y^2)\,s^2 +
        2s\,(r_{0,x}k_x + r_{0,y}k_y) + (r_{0,x}^2 + r_{0,y}^2 - \text{radius}^2) = 0, \f] whose
        smallest positive solution, if any, is the returned distance; kq2 is the leading
        coefficient \f$k_x^2+k_y^2\f$. */
    static double firstIntersectionCylinder(Vec r0, Vec k, double kq2, double radius)
    {
        if (std::abs(kq2) < EPS) return 0.;
        double b = r0.x() * k.x() + r0.y() * k.y();
        double c = square(r0.x()) + square(r0.y()) - square(radius);
        return smallestPositiveSolution(b / kq2, c / kq2);
    }

    /** This function determines the solutions of the intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the infinite cylinder described in firstIntersectionCylinder(), i.e.
        the two-root counterpart of that function: it solves the same defining equation and
        quadratic, but through distinctSolutions() rather than smallestPositiveSolution(), so that
        both roots are returned rather than only the smallest positive one; kq2 has the same
        meaning as there. The return value has the same meaning as for distinctSolutions(): 2 for
        two distinct solutions (stored in x1 and x2), or 0 if there are none, which includes the
        case where the ray is parallel to the z-axis (kq2 is zero). */
    static int distinctIntersectionsCylinder(Vec r0, Vec k, double kq2, double radius, double& x1, double& x2)
    {
        if (std::abs(kq2) < EPS) return 0;
        double b = r0.x() * k.x() + r0.y() * k.y();
        double c = square(r0.x()) + square(r0.y()) - square(radius);
        return distinctSolutions(b / kq2, c / kq2, x1, x2);
    }

    /** This function returns the distance to the intersection of the ray \f$({\bf{r}}_0,
        \hat{\bf{k}})\f$ with the horizontal plane, perpendicular to the z-axis, at the given z
        coordinate, or zero if the ray direction is parallel to the plane.

        Substituting \f${\bf{p}}(s) = {\bf{r}}_0 + s\hat{\bf{k}}\f$ into the plane's defining
        equation \f$p_z = z_\text{plane}\f$ and solving for \f$s\f$ gives \f[ s =
        \frac{z_\text{plane} - r_{0,z}}{k_z}. \f] */
    static double intersectionHorizontalPlane(Vec r0, Vec k, double zPlane)
    {
        if (std::abs(k.z()) < EPS) return 0.;
        return (zPlane - r0.z()) / k.z();
    }
};

////////////////////////////////////////////////////////////////////

#endif

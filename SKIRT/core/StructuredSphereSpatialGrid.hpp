/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STRUCTUREDSPHERESPATIALGRID_HPP
#define STRUCTUREDSPHERESPATIALGRID_HPP

#include "Mesh.hpp"
#include "SphereSpatialGrid.hpp"
class Random;

//////////////////////////////////////////////////////////////////////

/** StructuredSphereSpatialGrid is an abstract subclass of the SphereSpatialGrid class that
    represents a fully three-dimensional structured grid based on spherical coordinates. The grid
    is specified through three sets of grid points:

    - \f$N_r+1\f$ radial grid points \f$r_i, \,i=0,\ldots,N_r\f$, with \f$0\le r_\text{min} =
    r_0\f$, \f$r_i<r_{i+1}\f$, and \f$r_{N_r} = r_\text{max}\f$.

    - \f$N_\theta+1\f$ polar-inclination grid points \f$\theta_j, \,j=0,\ldots,N_\theta\f$ with
    \f$0=\theta_0\f$, \f$\theta_j<\theta_{j+1}\f$, and \f$\theta_{N_\theta}=\pi\f$.

    - \f$N_\varphi+1\f$ azimuthal grid points \f$\varphi_k, \,k=0,\ldots,N_\varphi\f$, with
    \f$-\pi=\varphi_0\f$, \f$0<\varphi_{k+1}-\varphi_k\le2\pi/3\f$, and
    \f$\varphi_{N_\varphi}=\pi\f$. The maximum limit on azimuth bin width is imposed to avoid
    confusion between the meridional half-planes when detecting cell border limits.

    \note The algorithm used by a path segment generator built on top of this class requires that
    the xy-plane \f$\theta=\pi/2\f$ is included in the polar grid. If this is not the case, this
    point is automatically added, increasing the number of polar bins by one.

    In total there are \f$N_{\text{cells}} = N_r\,N_\theta\,N_\varphi\f$ cells in the structured
    grid.

    This class provides the structured grid itself, including its geometric queries (cell volume,
    diagonal, containment, central and random positions, plotting functions), but leaves the path
    segment generator to concrete subclasses, since this depends heavily on whether anything
    besides the structured grid (such as embedded clumps) needs to be taken into account.
    Sphere3DSpatialGrid is a concrete subclass that adds nothing beyond the structured grid itself;
    ClumpySphericalSpatialGrid is a concrete subclass that superposes a set of spherical clumps on
    top of it. */
class StructuredSphereSpatialGrid : public SphereSpatialGrid
{
    ITEM_ABSTRACT(StructuredSphereSpatialGrid, SphereSpatialGrid,
                  "a structured 3D spatial grid in spherical coordinates")

        PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

        PROPERTY_ITEM(meshPolar, Mesh, "the bin distribution in the polar direction")
        ATTRIBUTE_DEFAULT_VALUE(meshPolar, "LinMesh")

        PROPERTY_ITEM(meshAzimuthal, Mesh, "the bin distribution in the azimuthal direction")
        ATTRIBUTE_DEFAULT_VALUE(meshAzimuthal, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up the structured grid. It precomputes and stores cosine values for the
        polar grid points and sine and cosine values for the azimuthal grid points, and it caches
        the nominal volume of each structured cell. It further ensures that the grid points conform
        to the requirements described in the class header. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 3. */
    int dimension() const override;

    /** This function returns the number of cells \f$N_r\,N_\theta\,N_\varphi\f$ in the structured
        grid. */
    int numCells() const override;

    /** This function returns the volume of the structured cell with index \f$m\f$, using the
        nominal volume cached during setup. The volume is calculated as \f[ V = \frac{1}{3}
        \left(r_{i+1}^3-r_i^3\right) \left(\cos\theta_j-\cos\theta_{j+1}\right)
        \left(\varphi_{k+1}-\varphi_k\right). \f]. */
    double volume(int m) const override;

    /** This function returns the "diagonal" of the structured cell with index \f$m\f$. It
        determines the 3D bin indices \f$i,j,k\f$ corresponding to the cell index \f$m\f$, and then
        calculates the distance between the outer/upper and inner/lower corners of the cell, i.e.
        between the points \f$\{ r_i,\theta_j,\varphi_k \}\f$ and \f$\{ r_{i+1},\theta_{j+1},
        \varphi_{k+1} \}\f$. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the structured cell that contains the position
        \f${\bf{r}}\f$, or -1 if the position is outside the grid. It determines the 3D bin indices
        \f$i,j,k\f$ of the cell containing the position and calculates the correct cell index based
        on these three numbers. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the structured cell with index \f$m\f$. It
        determines the 3D bin indices \f$i,j,k\f$ corresponding to the cell index \f$m\f$. The
        spherical coordinates of the central position are subsequently determined from \f[
        \begin{split} r &= \frac{r_i + r_{i+1}}{2} \\ \theta &= \frac{\theta_j + \theta_{j+1}}{2}
        \\ \varphi &= \frac{\varphi_k + \varphi_{k+1}}{2}. \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location in the structured cell with index \f$m\f$. It
        determines the 3D bin indices \f$i,j,k\f$ corresponding to the cell index \f$m\f$. Then a
        random radius \f$r\f$, a random inclination \f$\theta\f$, and a random azimuth
        \f$\varphi\f$ are determined using \f[ \begin{split} r &= \left( r_i^3 +
        {\cal{X}}_1\,(r_{i+1}^3-r_i^3) \right)^{1/3} \\ \cos\theta &= \cos\theta_j + {\cal{X}}_2\,
        (\cos\theta_{j+1}-\cos\theta_j) \\ \varphi &= \varphi_k + {\cal{X}}_3\,
        (\varphi_{k+1}-\varphi_k), \end{split} \f] with \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and
        \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

protected:
    /** This function writes the intersection of the structured grid with the xy plane to the
        specified SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the structured grid with the xz plane to the
        specified SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the structured grid with the yz plane to the
        specified SpatialGridPlotFile object. */
    void write_yz(SpatialGridPlotFile* outfile) const override;

    /** This function writes 3D information for the structured grid to the specified
        SpatialGridPlotFile object. A subclass that superposes additional geometry on the
        structured grid may override this function to overlay that geometry, for example by
        calling this implementation first and then adding its own content. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

    /** This function returns the relative-to-maxRadius scale factor used, together with the mesh's
        own natural scale (the width of the innermost bin or, if there is a hole, the hole's
        radius), to determine the small epsilon used for progressing a path across a structured
        cell boundary. The default implementation returns a value that has proven robust for a
        plain structured grid. A subclass that superposes additional internal geometry (such as
        embedded clumps) on the structured grid may need to override this with a more conservative
        (larger) value, since paths in that subclass cross more boundaries of more varied kinds
        within a given region. */
    virtual double meshEpsilonScale() const;

    /** This function returns the cell index corresponding to the radial index \f$i\f$, the polar
        index \f$j\f$, and azimuthal index \f$k\f$. The correspondence is
        \f$m=k+j\,N_\varphi+i\,N_\theta\,N_\varphi\f$. */
    int index(int i, int j, int k) const;

    /** This function writes the axisymmetric structured-grid lines (spheres and cones) intersected
        by any meridional plane to the specified SpatialGridPlotFile object. Because the structured
        grid is itself axisymmetric, this single function serves both the xz and yz views; a
        subclass overriding write_xz() or write_yz() to overlay additional, non-axisymmetric
        content can call this function to obtain the structured grid's own contribution. */
    void writeMeridionalStructure(SpatialGridPlotFile* outfile) const;

    /** This function obtains the spherical coordinates for the corners of the structured cell with
        index \f$m\f$. It determines the radial, polar and azimuthal bin indices \f$i\f$, \f$j\f$
        and \f$k\f$ corresponding to the cell index \f$m\f$ using the formulae \f[ \begin{split} i
        &= \lfloor m/(N_\theta\,N_\varphi) \rfloor \\ j &= \lfloor m/N_\varphi \rfloor
        \,\text{mod}\, N_\theta \\ k &= m\,\text{mod}\,N_\varphi. \end{split} \f]

        If all of the resulting bin indices are within range, the function stores the corresponding
        cell corner coordinates in the provided arguments and returns true. If any of the indices
        are out of range, the function returns false and the contents of the provided arguments
        remains unchanged. */
    bool getCoords(int m, double& rmin, double& thetamin, double& phimin, double& rmax, double& thetamax,
                   double& phimax) const;

    //======================== Data Members ========================

protected:
    int _Nr{0};
    int _Ntheta{0};
    int _Nphi{0};
    int _Ncells{0};
    Array _rv;
    Array _thetav;
    Array _phiv;
    Array _cv;
    Array _sinv;
    Array _cosv;
    double _eps{0.};

    // the nominal volume of each structured cell (index m, 0 <= m < _Ncells), cached during setup;
    // a subclass that reduces some cells' volume by overlapping internal geometry may adjust this
    // array in place after calling this class's setupSelfAfter()
    vector<double> _cellVolume;
};

//////////////////////////////////////////////////////////////////////

#endif

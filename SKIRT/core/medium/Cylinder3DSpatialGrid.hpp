/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CYLINDER3DSPATIALGRID_HPP
#define CYLINDER3DSPATIALGRID_HPP

#include "Array.hpp"
#include "CylinderSpatialGrid.hpp"
#include "Mesh.hpp"
class CylindricalCell;

////////////////////////////////////////////////////////////////////

/** The Cylinder3DSpatialGrid class is subclass of the CylinderSpatialGrid class, and represents a
    fully three-dimensional spatial grid based on cylindrical coordinates. The grid is specified
    through three sets of grid points:

    - \f$N_R+1\f$ radial grid points \f$R_i, \,i=0,\ldots,N_R\f$, with \f$0\le R_\text{min} =
    R_0\f$, \f$R_i<R_{i+1}\f$, and \f$R_{N_R} = R_\text{max}\f$.

    - \f$N_\varphi+1\f$ azimuthal grid points \f$\varphi_j, \,j=0,\ldots,N_\varphi\f$, with
    \f$-\pi=\varphi_0\f$, \f$0<\varphi_{j+1}-\varphi_j\le2\pi/3\f$, and
    \f$\varphi_{N_\varphi}=\pi\f$. The maximum limit on azimuth bin width is imposed to avoid
    confusion between the meridional half-planes when detecting cell border limits.

    - \f$N_z+1\f$ vertical grid points \f$z_k, \,k=0,\ldots,N_z\f$, with \f$z_\text{min} = z_0\f$,
    \f$z_i<z_{i+1}\f$, and \f$z_{N_z} = z_\text{max}\f$.

    In total there are \f$N_{\text{cells}} = N_R\,N_\varphi\,N_z\f$ cells in the grid. */
class Cylinder3DSpatialGrid : public CylinderSpatialGrid
{
    ITEM_CONCRETE(Cylinder3DSpatialGrid, CylinderSpatialGrid, "a 3D spatial grid in cylindrical coordinates")

        PROPERTY_ITEM(meshRadial, Mesh, "the bin distribution in the radial direction")
        ATTRIBUTE_DEFAULT_VALUE(meshRadial, "LinMesh")

        PROPERTY_ITEM(meshAzimuthal, Mesh, "the bin distribution in the azimuthal direction")
        ATTRIBUTE_DEFAULT_VALUE(meshAzimuthal, "LinMesh")

        PROPERTY_ITEM(meshZ, Mesh, "the bin distribution in the vertical direction")
        ATTRIBUTE_DEFAULT_VALUE(meshZ, "LinMesh")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function sets up a number of data members that depend on the Mesh objects configured
        for this grid. It precomputes and stores sine and cosine values for the azimuthal grid
        points, and verifies that the azimuth bins are not too wide. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 3. */
    int dimension() const override;

    /** This function returns the number of cells \f$N_R\,N_\varphi\,N_z\f$ in the grid. */
    int numCells() const override;

    /** This function returns the volume of the cell with index \f$m\f$. It determines the 3D bin
        indices \f$i,j,k\f$ corresponding to the index \f$m\f$, and then calculates the volume as
        \f$\frac{1}{2} (R_{i+1}^2-R_i^2) (\varphi_{j+1}-\varphi_j) (z_{k+1}-z_k)\f$. */
    double volume(int m) const override;

    /** This function returns the "diagonal" of the cell with index \f$m\f$. It determines the 3D
        bin indices \f$i,j,k\f$ corresponding to the index \f$m\f$, and then calculates the
        distance between the outer/upper and inner/lower corners of the cell, i.e. between the
        points \f$\{ R_i,\varphi_j,z_k \}\f$ and \f$\{ R_{i+1},\varphi_{j+1},z_{k+1} \}\f$. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. It determines the indices \f$i,j,k\f$ of the cell containing the position
        and calculates the correct cell index based on these numbers. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. It determines
        the 3D bin indices \f$i,j,k\f$ corresponding to the index \f$m\f$, and then calculates the
        cylindrical coordinates of the central position using \f[ \begin{split} R &= \frac{R_i +
        R_{i+1}}{2} \\ \varphi &= \frac{\varphi_j + \varphi_{j+1}}{2} \\ z &= \frac{z_k +
        z_{k+1}}{2} \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random position in the cell with index \f$m\f$. It determines the
        3D bin indices \f$i,j,k\f$ corresponding to the index \f$m\f$, and then generates the
        cylindrical coordinates of a random position using \f[ \begin{split} R &= \sqrt{R_i^2 +
        {\cal{X}}_1\,(R_{i+1}^2-R_i^2)} \\ \varphi &= \varphi_j + {\cal{X}}_2\,
        (\varphi_{j+1}-\varphi_j) \\ z &= z_k + {\cal{X}}_3\, (z_{k+1}-z_k), \end{split} \f] with
        \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for a 3D cylindrical grid, implemented as a
        private PathSegmentGenerator subclass.

        The ray equation of the path can be written as \f[\begin{cases} x = r_\text{x} +
        k_\text{x}s \\ y = r_\text{y} + k_\text{y}s \\ z = r_\text{z} + k_\text{z}s \\ \end{cases}
        \quad \text{with} \;s>0.\f]

        Intersection with a vertical cylinder centered on the origin with equation \f$x^2 +y^2 =
        R_*^2\f$ yields a quadratic equation of the form \f$s^2+2bs+c=0\f$ with \f[\begin{aligned}
        b &= \frac{r_\text{x}k_\text{x} + r_\text{y}k_\text{y}} {k_\text{x}^2 + k_\text{y}^2} \\ c
        &= \frac{r_\text{x}^2 + r_\text{y}^2 - R_*^2} {k_\text{x}^2 + k_\text{y}^2} \\
        \end{aligned}\f]

        Intersection with a meridional plane with equation \f$\sin\varphi_* x = \cos\varphi_* y\f$
        yields \f[ s = -\;\frac{r_\text{x}\sin\varphi_* - r_\text{y}\cos\varphi_*}
        {k_\text{x}\sin\varphi_* - k_\text{y}\cos\varphi_*} \f] The requirement that
        \f$\varphi_{j+1}-\varphi_j\le2\pi/3\f$ ensures that the intersection point with the
        reflected half-plane is always more distant than the other cell boundaries thus that
        phantom point is automatically ignored.

        Intersection with a horizontal plane with equation \f$z=z_*\f$ easily yields \f$s =
        (z_*-r_\text{z})/k_\text{z}\f$.

        The segment generator progresses the starting point of the path through the grid along the
        path's direction. For each step along the way, it calculates the distances to the
        intersections with all candidate borders of the current cell, and then selects the nearest
        intersection point. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the xz plane to the specified
        SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the yz plane to the specified
        SpatialGridPlotFile object. */
    void write_yz(SpatialGridPlotFile* outfile) const override;

    /** This function writes 3D information for all or part of the cells in the grid structure to
        the specified SpatialGridPlotFile object. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

private:
    /** This function returns the index \f$m\f$ corresponding to the 3D bin indices \f$i\f$,
        \f$j\f$ and \f$k\f$. The correspondence is \f$m=k+j\,N_z+i\,N_\varphi\,N_z\f$. */
    int index(int i, int j, int k) const;

    /** This function obtains the cylindrical coordinates for the corners of the cell with index
        \f$m\f$. It determines the 3D bin indices \f$i,j,k\f$ corresponding to the index \f$m\f$
        using the formulae \f[ \begin{split} i &= \lfloor m/(N_\varphi\,N_z) \rfloor \\ j &=
        \lfloor m/N_z \rfloor \,\text{mod}\, N_\varphi \\ k &= m\,\text{mod}\,N_z. \end{split} \f]

        If all of the resulting bin indices are within range, the function stores the corresponding
        cell corner coordinates in the provided arguments and returns true. If any of the indices
        are out of range, the function returns false and the contents of the provided arguments
        remains unchanged. */
    bool getCoords(int m, double& Rmin, double& phimin, double& zmin, double& Rmax, double& phimax, double& zmax) const;

    //======================== Data Members ========================

private:
    int _NR{0};
    int _Nphi{0};
    int _Nz{0};
    int _Ncells{0};
    Array _Rv;
    Array _phiv;
    Array _zv;
    Array _sinv;
    Array _cosv;
    double _eps{0.};
    bool _hasHole{false};

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif

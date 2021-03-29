/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIMESHMEDIUM_HPP
#define VORONOIMESHMEDIUM_HPP

#include "MeshMedium.hpp"
#include "VoronoiMeshInterface.hpp"

////////////////////////////////////////////////////////////////////

/** A VoronoiMeshMedium instance represents a transfer medium with a spatial density distribution
    described by a list of sites generating a Voronoi tesselation of a cubodail domain. The data is
    usually extracted from a cosmological simulation snapshot, and it must be provided in a column
    text file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this VoronoiMeshMedium instance:

    \f[ x\,(\mathrm{pc}) \quad y\,(\mathrm{pc}) \quad z\,(\mathrm{pc}) \quad \{\,
    \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \quad [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) ] \quad [
    B_x\,(\mu\mathrm{G}) \quad B_y\,(\mu\mathrm{G}) \quad B_z\,(\mu\mathrm{G}) ] \quad [
    \dots\text{mix family params}\dots ] \f]

    The first three columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of the Voronoi site
    (i.e. the location defining a particular Voronoi cell).

    Depending on the value of the \em massType option, the fourth column lists the average mass
    density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or the
    integrated number density \f$N\f$ for the cell corresponding to the site. This quantity is
    multiplied by the value of the \em massFraction option.

    If the \em importMetallicity option is enabled, the next column specifies a "metallicity"
    fraction, which is multiplied with the mass/density column to obtain the actual mass/density of
    the cell corresponding to the site. If the \em importTemperature option is enabled, the next
    column specifies a temperature. If this temperature is higher than the value of the \em
    maxTemperature option, the mass and density for the site are set to zero, regardless of the
    mass or density specified in the fourth column. If the \em importTemperature option is
    disabled, or the maximum temperature value is set to zero, such a cutoff is not applied.

    If both the \em importMetallicity and \em importTemperature options are enabled, this leads to
    the following expression for the density of an imported site (or a simular formula for the
    other mass quantity types):

    \f[ \rho_\mathrm{imported} = \begin{cases} f_\mathrm{mass}\,Z\,\rho & \mathrm{if}\;
    T<T_\mathrm{max} \;\mathrm{or}\; T_\mathrm{max}=0 \\ 0 & \mathrm{otherwise} \end{cases} \f]

    If the \em importVelocity option is enabled, the subsequent three columns specify the
    \f$v_x\f$, \f$v_y\f$, \f$v_z\f$ velocity components of the site (considered as the bulk
    velocity for the mass represented by the site).

    If the \em importMagneticField option is enabled, the subsequent three columns specify the
    \f$B_x\f$, \f$B_y\f$, \f$B_z\f$ magnetic field vector components for the cell.

    Finally, if the \em importVariableMixParams option is enabled, the remaining columns specify
    the parameters used by the configured material mix family to select a particular material mix
    for the cell.

    <b>Avoiding construction of the Voronoi tessellation</b>

    The algorithm used by the VoronoiMeshSnapshot class for constructing Voronoi tessellations
    sometimes fails (for example, when generating sites are too close to each other). In those
    cases, it can be desirable to forego the construction of the Voronoi tessellation and instead
    use a nearest neighbor search for sampling the density distribution. This is possible as long
    as the full tessellation is not needed for other purposes, such as to perform radiative
    transfer or to generate random positions drawn from the density distribution. An important use
    case is when the medium density distribution defined on the Voronoi grid is resampled to an
    octree grid to perform radiative transfer. However, the octree subdivision algorithm requires
    the total mass of the medium in addition to the density samples. Without the full Voronoi
    tessellation, it is impossible to calculate the cell volume that would allow conversion between
    cell density and mass.

    To enable this use case, the \em massType option can be set to include both mass density and
    volume-integrated mass (or both number density and volume-integrated number) in the imported
    data file. Thus, if the \em massType option is set to one of these choices, the import file
    must include two columns (mass density \f$\rho\f$ plus integrated mass \f$M\f$, or number
    density \f$n\f$ plus integrated number density \f$N\f$ ). Furthermore, if the simulation
    configuration allows it (i.e. the full Voronoi tessellation is not needed for other purposes),
    the VoronoiMeshSnapshot class will use the information in these two columns to calculate the
    cell volumes and the total medium mass, and it will forego construction of the Voronoi
    tessellation. */
class VoronoiMeshMedium : public MeshMedium, public VoronoiMeshInterface
{
    ITEM_CONCRETE(VoronoiMeshMedium, MeshMedium, "a transfer medium imported from data represented on a Voronoi mesh")
        ATTRIBUTE_TYPE_INSERT(VoronoiMeshMedium, "VoronoiMeshInterface")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new VoronoiMeshSnapshot object, calls its open() function,
        passes it the domain extent configured by the user, configures it to import a mass or a
        density column, and finally returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;

    //=================== Other functions ==================

protected:
    /** This function implements the VoronoiMeshInterface interface. It returns a pointer to the
        Voronoi mesh snapshot maintained by this geometry. */
    VoronoiMeshSnapshot* voronoiMesh() const override;

    //===================== Data members ====================

private:
    // an extra pointer to our snapshot used to implement VoronoiMeshInterface (ownership is passed to base class)
    VoronoiMeshSnapshot* _voronoiMeshSnapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif

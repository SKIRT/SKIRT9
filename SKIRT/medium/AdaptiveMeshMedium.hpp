/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHMEDIUM_HPP
#define ADAPTIVEMESHMEDIUM_HPP

#include "AdaptiveMeshInterface.hpp"
#include "MeshMedium.hpp"

////////////////////////////////////////////////////////////////////

/** An AdaptiveMeshMedium instance represents a transfer medium with a spatial density distribution
    described by an Adaptive Mesh Refinement (AMR) grid partitioning a cuboidal domain. The data is
    usually extracted from a cosmological simulation snapshot, and it must be provided in a column
    text file formatted as described below.

    Refer to the description of the AdaptiveMeshSnapshot class for information on the structure of
    an adaptive mesh and on how to represent it in text column file format. Refer to the
    description of the TextInFile class for information on overall formatting and on how to include
    header lines specifying the units for each column in the input file. In case the input file has
    no unit specifications, the default units mentioned below are used instead. The number of
    columns expected in the input file depends on the options configured by the user for this
    AdaptiveMeshMedium instance:

    \f[ \{\, \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \quad [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) ] \quad [
    B_x\,(\mu\mathrm{G}) \quad B_y\,(\mu\mathrm{G}) \quad B_z\,(\mu\mathrm{G}) ] \quad [
    \dots\text{mix family params}\dots ] \f]

    Depending on the value of the \em massType option, the first column lists the average mass
    density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or the
    integrated number density \f$N\f$ for the cell. This quantity is multiplied by the value of the
    \em massFraction option.

    If the \em importMetallicity option is enabled, the next column specifies a "metallicity"
    fraction, which is multiplied with the mass/density column to obtain the actual mass/density of
    the cell. If the \em importTemperature option is enabled, the next column specifies a
    temperature. If this temperature is higher than the value of the \em maxTemperature option, the
    mass and density for the cell are set to zero, regardless of the mass or density specified in
    the first column. If the \em importTemperature option is disabled, or the maximum temperature
    value is set to zero, such a cutoff is not applied.

    If both the \em importMetallicity and \em importTemperature options are enabled, this leads to
    the following expression for the density of an imported cell (or a simular formula for the
    other mass quantity types):

    \f[ \rho_\mathrm{imported} = \begin{cases} f_\mathrm{mass}\,Z\,\rho & \mathrm{if}\;
    T<T_\mathrm{max} \;\mathrm{or}\; T_\mathrm{max}=0 \\ 0 & \mathrm{otherwise} \end{cases} \f]

    If the \em importVelocity option is enabled, the subsequent three columns specify the
    \f$v_x\f$, \f$v_y\f$, \f$v_z\f$ velocity components of the cell (considered as the average bulk
    velocity for the mass in the cell).

    If the \em importMagneticField option is enabled, the subsequent three columns specify the
    \f$B_x\f$, \f$B_y\f$, \f$B_z\f$ magnetic field vector components for the cell.

    Finally, if the \em importVariableMixParams option is enabled, the remaining columns specify
    the parameters used by the configured material mix family to select a particular material mix
    for the cell. */
class AdaptiveMeshMedium : public MeshMedium, public AdaptiveMeshInterface
{
    ITEM_CONCRETE(AdaptiveMeshMedium, MeshMedium,
                  "a transfer medium imported from data represented on an adaptive mesh (AMR grid)")
        ATTRIBUTE_TYPE_INSERT(AdaptiveMeshMedium, "AdaptiveMeshInterface")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new AdaptiveMeshSnapshot object, calls its open() function,
        passes it the domain extent configured by the user, configures it to import a mass or a
        density column, and finally returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;

    //=================== Other functions ==================

protected:
    /** This function implements the AdaptiveMeshInterface interface. It returns a pointer to the
        adaptive mesh snapshot maintained by this geometry. */
    AdaptiveMeshSnapshot* adaptiveMesh() const override;

    //===================== Data members ====================

private:
    // an extra pointer to our snapshot used to implement AdaptiveMeshInterface (ownership is passed to base class)
    AdaptiveMeshSnapshot* _adaptiveMeshSnapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif

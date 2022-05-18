/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARTICLEMEDIUM_HPP
#define PARTICLEMEDIUM_HPP

#include "ImportedMedium.hpp"
#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** A ParticleMedium instance represents a transfer medium with a spatial density distribution
    (and, optionally, other properties) described by a list of smoothed particles. The particle
    data is usually extracted from a cosmological simulation snapshot, and it must be provided in a
    column text file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this ParticleMedium instance:

    \f[ x\,(\mathrm{pc}) \quad y\,(\mathrm{pc}) \quad z\,(\mathrm{pc}) \quad h\,(\mathrm{pc}) \quad
    \{\, M\,(\mathrm{M}_\odot) \;\;|\;\; N\,(1) \,\} \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \quad [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) ] \quad [
    B_x\,(\mu\mathrm{G}) \quad B_y\,(\mu\mathrm{G}) \quad B_z\,(\mu\mathrm{G}) ] \quad [
    \dots\text{mix family params}\dots ] \f]

    The first three columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of the particle, the
    fourth column is the particle smoothing length \f$h\f$.

    The fifth column is the volume-integrated mass \f$M\f$ or number \f$N\f$ of the particle,
    depending on the value of the \em massType option. This mass or number is multiplied by the
    value of the \em massFraction option. If the \em importMetallicity option is enabled, the next
    column specifies a "metallicity" fraction, which is multiplied with the mass or number column
    to obtain the actual mass of the particle. If the \em importTemperature option is enabled, the
    next column specifies a temperature. If this temperature is higher than the value of the \em
    maxTemperature option, the particle is ignored. If the \em importTemperature option is
    disabled, or the maximum temperature value is set to zero, the particle is never ignored.

    If both the \em importMetallicity and \em importTemperature options are enabled, this leads to
    the following expression for the mass of an imported particle:

    \f[ M_\mathrm{imported} = \begin{cases} f_\mathrm{mass}\,Z\,M & \mathrm{if}\; T<T_\mathrm{max}
    \;\mathrm{or}\; T_\mathrm{max}=0 \\ 0 & \mathrm{otherwise} \end{cases} \f]

    If the \em importVelocity option is enabled, the subsequent three columns specify the
    \f$v_x\f$, \f$v_y\f$, \f$v_z\f$ velocity components of the particle (considered as the bulk
    velocity for the mass represented by the particle).

    If the \em importMagneticField option is enabled, the subsequent three columns specify the
    \f$B_x\f$, \f$B_y\f$, \f$B_z\f$ magnetic field vector components for the cell.

    Finally, if the \em importVariableMixParams option is enabled, the remaining columns specify
    the parameters used by the configured material mix family to select a particular material mix
    for the particle. */
class ParticleMedium : public ImportedMedium
{
    /** The enumeration type indicating the type of mass quantity to be imported. */
    ENUM_DEF(MassType, Mass, Number)
        ENUM_VAL(MassType, Mass, "mass (volume-integrated)")
        ENUM_VAL(MassType, Number, "number (volume-integrated)")
    ENUM_END()

    ITEM_CONCRETE(ParticleMedium, ImportedMedium, "a transfer medium imported from smoothed particle data")

        PROPERTY_ITEM(smoothingKernel, SmoothingKernel, "the kernel for interpolating the smoothed particles")
        ATTRIBUTE_DEFAULT_VALUE(smoothingKernel, "CubicSplineSmoothingKernel")
        ATTRIBUTE_DISPLAYED_IF(smoothingKernel, "Level2")

        PROPERTY_ENUM(massType, MassType, "the type of mass quantity to be imported")
        ATTRIBUTE_DEFAULT_VALUE(massType, "Mass")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new ParticleSnapshot object, calls its open() function,
        configures it to import a mass column, passes the smoothing kernel selected by the user to
        it, and finally returns a pointer to the object. Ownership of the Snapshot object is
        transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif

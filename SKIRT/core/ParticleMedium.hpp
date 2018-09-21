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
    input file should contain 5 to 10 columns, depending on the options configured by the user for
    this ParticleMedium instance:

    \f[ x\,(\mathrm{pc}) \quad y\,(\mathrm{pc}) \quad z\,(\mathrm{pc}) \quad h\,(\mathrm{pc}) \quad
    M\,(\mathrm{M}_\odot) \quad [Z\,(1)] \quad [T\,(\mathrm{K})] \quad [ v_x\,(\mathrm{km/s}) \quad
    v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) ] \f]

    The first three columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of the particle, the
    fourth column is the particle smoothing length \f$h\f$.

    The fifth column is the mass \f$M\f$ of the particle, which is multiplied by the value of the
    \em massFraction option. If the \em importMetallicity option is enabled, the next column
    specifies a "metallicity" fraction, which is multiplied with the mass column to obtain the
    actual mass of the particle. If the \em importTemperature option is enabled, the next column
    specifies a temperature. If this temperature is higher than the value of the \em maxTemperature
    option, the particle is ignored. If the \em importTemperature option is disabled, or the
    maximum temperature value is set to zero, the particle is never ignored.

    If both the \em importMetallicity and \em importTemperature options are enabled, this leads to
    the following expression for the mass of an imported particle:

    \f[ M_\mathrm{imported} = \begin{cases} f_\mathrm{mass}\,Z\,M & \mathrm{if}\; T<T_\mathrm{max}
    \;\mathrm{or}\; T_\mathrm{max}=0 \\ 0 & \mathrm{otherwise} \end{cases} \f]

    If the \em importVelocity option is enabled, the final three columns specify the \f$v_x\f$,
    \f$v_y\f$, \f$v_z\f$ velocity components of the particle (considered as the bulk velocity for
    the mass represented by the particle). */
class ParticleMedium : public ImportedMedium
{
    ITEM_CONCRETE(ParticleMedium, ImportedMedium, "a transfer medium imported from smoothed particle data")

    PROPERTY_ITEM(smoothingKernel, SmoothingKernel, "the kernel for interpolating the smoothed particles")
        ATTRIBUTE_DEFAULT_VALUE(smoothingKernel, "CubicSplineSmoothingKernel")
        ATTRIBUTE_DISPLAYED_IF(smoothingKernel, "Level2")

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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PLANARCUSTOMSTATECUTSPROBE_HPP
#define PLANARCUSTOMSTATECUTSPROBE_HPP

#include "AbstractPlanarCutsProbe.hpp"

////////////////////////////////////////////////////////////////////

/** For each medium in the simulation that has one or more custom medium state variables (as
    requested by the associated material mix), PlanarCustomStateCutsProbe outputs a FITS file with
    cuts through the values of those custom state variables along three planes parallel to the
    coordinate planes. The offset of each cut plane from the corresponding coordinate plane can be
    configured by the user (and is zero by default). The field of view of each cut covers the
    extent of the spatial grid in the simulation in the relevant directions. The number of pixels
    in each direction can be configured by the user as well.

    The files are named <tt>prefix_probe_customstate_XX_N.fits</tt> where XX indicates the
    orientation of the cut and N is replaced with the zero-based index of the medium in the
    configuration (i.e. in the ski file). Each file contains an image frame for each of the custom
    state variables of the corresponding medium, in the order in which those variables are
    requested by the associated material mix.

    \note The current implementation assumes that all custom state variables for a given medium
    component are the same physical quantity type and thus also have the same output units. In
    principle this restriction could be lifted but in that case it is unclear where to put the unit
    information in the FITS header. */
class PlanarCustomStateCutsProbe : public AbstractPlanarCutsProbe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
    ENUM_END()

    ITEM_CONCRETE(PlanarCustomStateCutsProbe, AbstractPlanarCutsProbe,
                  "cuts of the custom medium state along planes parallel to the coordinate planes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(PlanarCustomStateCutsProbe, "Level2&CustomMediumState&SpatialGrid")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "when to probe the medium state")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "HasDynamicState")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function outputs FITS files with custom medium state cuts in a plane parallel to the
        coordinate plane indicated by the boolean "direction" arguments \em xd, \em yd, and \em zd,
        exactly two of which must be true. The arguments \em xc, \em yc, and \em zc specify the
        position of the cuts, and the arguments \em Nx, \em Ny, and \em Nz specify the number of
        pixels in each direction. */
    static void writeCustomStateCuts(Probe* probe, bool xd, bool yd, bool zd, double xc, double yc, double zc, int Nx,
                                     int Ny, int Nz);

    /** This function performs probing after setup. It produces output only if the \em
        probeAfter property is set to Setup. */
    void probeSetup() override;

    /** This function performs probing after all photon packets have been emitted and detected. It
        produces output only if the \em probeAfter property is set to Run. */
    void probeRun() override;

private:
    /** This function performs the probing; it is called from probeSetup() or probeRun() depending
        on the value of the \em probeAfter property. */
    void probe();
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTGRAINPOPULATIONSPROBE_HPP
#define DUSTGRAINPOPULATIONSPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** DustGrainPopulationsProbe outputs column text files listing information on the dust grain
    populations of the applicable dust media configured in the simulation. For each medium
    component, the probe retrieves a representative material mix (the mix at the origin of the
    model coordinate system). If the material mix offers MultiGrainDustMix capabilities (i.e. it is
    described by one or more individual dust grain populations), the probe creates a file with
    information on each of the grain populations in that mix. The files are named
    <tt>prefix_probe_grainpops_N.fits</tt> where N is replaced with the zero-based index of the
    medium in the configuration (i.e. in the ski file).

    The columns in each text file list the following properties for each grain population: the
    zero-based population index; a short name identifying the type of grain material in the
    population; the dust mass in the population as a percentage of the total, as a mass per
    hydrogen atom, and as a mass per hydrogen mass; and the size range of the grains in the
    population. */
class DustGrainPopulationsProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(DustGrainPopulationsProbe, SpecialtyProbe,
                  "properties: dust grain population mass and size information")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DustGrainPopulationsProbe, "Level2&Medium&MultiGrainDustMix")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

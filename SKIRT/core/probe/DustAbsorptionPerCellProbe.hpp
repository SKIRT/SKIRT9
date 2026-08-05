/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTABSORPTIONPERCELLPROBE_HPP
#define DUSTABSORPTIONPERCELLPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** DustAbsorptionPerCellProbe outputs a column text file (named <tt>prefix_probe_Labs.dat</tt>)
    listing the spectral luminosity \f$L_\lambda^\text{abs}\f$ absorbed by dust for each cell in the
    spatial grid of the simulation. It is calculated from the mean intensity \f$J_\lambda\f$ as
    follows,

    \f[ (L_\lambda^\text{abs})_{\ell,m} = (J_\lambda)_{\ell,m} \,4\pi\,V_m \sum_h^\text{dust}
    \kappa_{\ell,h}^{abs} \rho_{m,h}\f]

    where \f$\ell\f$ is the index of the wavelength bin, \f$m\f$ is the spatial cell index,
    \f$V_m\f$ is the volume of the spatial cell, and the sum represents the dust opacity in the
    cell.

    Specifically, the output file contains a line for each cell in the spatial grid of the
    simulation. The first columns specifies the cell index, and subsequent columns list the
    spectral luminosity absorbed by dust for each bin in the wavelength grid returned by the
    Configuration::radiationFieldWLG() function.

    The probe offers an option to output a separate text column file with details on the radiation
    field wavelength grid. For each wavelength bin, the file lists the characteristic wavelength,
    the wavelength bin width, and the left and right borders of the bin. */
class DustAbsorptionPerCellProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(DustAbsorptionPerCellProbe, SpecialtyProbe,
                  "specialty: spectral luminosity absorbed by dust for each spatial cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DustAbsorptionPerCellProbe, "Level3&DustMix&SpatialGrid&RadiationField")

        PROPERTY_BOOL(writeWavelengthGrid, "output a text file with the radiation field wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(writeWavelengthGrid, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns the enumeration \c Run indicating that probing for this probe should
        be performed at the end of the simulation. */
    When when() const override;

    //======================== Other Functions =======================

protected:
    /** This function performs probing after all photon packets have been emitted and detected. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

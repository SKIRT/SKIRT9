/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIVITYPROBE_HPP
#define DUSTEMISSIVITYPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** DustEmissivityProbe outputs column text files tabulating the emissivity spectra for the
    representative dust mix of each medium in the medium system (i.e. the dust mix at the model
    coordinate origin), using the dust emissivity calculator configured for the simulation, and
    assuming the dust would be embedded in a range of predefined input fields.

    The output emissivity spectra are discretized on the emission spectrum wavelength grid as
    returned by the Configuration::dustEmissionWLG() function. The input radiation fields used
    for the calculations are discretized on the radiation field wavelength grid as returned by the
    Configuration::radiationFieldWLG() function. As a result, this probe requires the configuration
    of the simulation to include dust emission, so that a dust emissivity calculator, a radiation
    field wavelength grid, and an emission spectrum wavelength grid have been properly configured.

    The built-in embedding fields include scaled versions of the solar neighborhood interstellar
    radiation field presented by Mathis et al. (1983, A&A, 128, 212), which is essentially a sum of
    three diluted blackbodies with a UV extension, and a set of diluted blackbodies at varying
    temperatures.

    A separate output file is written for each embedding field. Each file contains a line per
    wavelength. The first column lists the wavelength, amd the subsequent columns list the
    emmissivity at that wavelength for each of the medium components that actually contain dust.
    The emmissivity values are always given in units of Watt per steradian per hydrogen atom
    ("W/sr/H"), regardless of the configured output unit style.

    The probe offers an option to output a separate text column file with additional details on the
    emission spectrum wavelength grid. For each wavelength bin, the file lists the characteristic
    wavelength, the wavelength bin width, and the left and right borders of the bin. */
class DustEmissivityProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(DustEmissivityProbe, SpecialtyProbe,
                  "properties: emissivity spectrum for each dust mix in a range of standard fields")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DustEmissivityProbe, "Level3&DustMix&DustEmission")

        PROPERTY_BOOL(writeWavelengthGrid, "output a text file with the emission spectrum wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(writeWavelengthGrid, "false")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

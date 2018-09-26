/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPTICALPROPERTIESPROBE_HPP
#define OPTICALPROPERTIESPROBE_HPP

#include "Probe.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

/** OpticalPropertiesProbe outputs column text files listing the key optical properties for the
    media configured in the simulation, discretized on a specified wavelength grid or on the
    default instrument wavelength grid. For each medium component, the probe retrieves a
    representative material mix (the mix at the origin of the model coordinate system) and creates
    a file with the key optical properties for that mix. The files are named
    <tt>prefix_probe_opticalprops_N.fits</tt> where N is replaced with the zero-based index of the
    medium in the configuration (i.e. in the ski file).

    The first line in each file indicates the fundamental material type (i.e. dust, electrons, or
    gas) and lists the mass per material entity (hydrogen atom or electron) \f$\mu\f$, which can be
    used to convert the cross sections to mass coefficients.

    The columns list the key aggregate properties for the material mix as a function of wavelength.
    In order from left ot right, the wavelength \f$\lambda\f$; the total extinction, absorption,
    and scattering cross sections per hydrogen atom (or electron) \f$\varsigma_\lambda^\text{ext},
    \varsigma_\lambda^\text{abs}, \varsigma_\lambda^\text{sca}\f$; the total extinction,
    absorption, and scattering mass coefficients \f$\kappa_\lambda^\text{ext},
    \kappa_\lambda^\text{abs}, \kappa_\lambda^\text{sca}\f$; the corresponding scattering albedo
    \f$\varpi_\lambda\f$; and the mean asymmetry parameter \f$g_\lambda\f$ (or zero if not
    available). See the MaterialMix class for more information. */
class OpticalMaterialPropertiesProbe : public Probe
{
    ITEM_CONCRETE(OpticalMaterialPropertiesProbe, Probe, "aggregate optical material properties for each medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(OpticalMaterialPropertiesProbe, "Medium")

    PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for the optical properties")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultInstrumentWavelengthGrid")
        ATTRIBUTE_DISPLAYED_IF(wavelengthGrid, "Level2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after setup. */
    void probeSetup() override;
};

////////////////////////////////////////////////////////////////////

#endif

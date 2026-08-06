
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIXFAMILY_HPP
#define XRAYIONICGASMIXFAMILY_HPP

#include "MaterialMixFamily.hpp"
#include "XRayIonicGasMix.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the XRayIonicGasMixFamily class represents a family of dust mixes that is
    specified as part of the configuration. Specifically, a property of this class holds a
    user-configurable list of dust mixes representing the family. The family requires a single
    parameter value to select a family member, corresponding to the zero-based index in the
    configured list of dust mixes. The floating point parameter value is rounded to the nearest
    integer and subsequently clipped to be in range. */
class XRayIonicGasMixFamily : public MaterialMixFamily
{
    ENUM_DEF(ElectronScattering, None, Free, FreeWithPolarization)
        ENUM_VAL(ElectronScattering, None, "ignore electron")
        ENUM_VAL(ElectronScattering, Free, "use free-electron Compton scattering for all electrons")
        ENUM_VAL(ElectronScattering, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
    ENUM_END()

    ITEM_CONCRETE(XRayIonicGasMixFamily, MaterialMixFamily, "a family of ionic mixes for each cell")

        PROPERTY_STRING(ions, "the names of the ions for each element (e.g. H,He+,Li+1,..)")

        PROPERTY_ENUM(electronScattering, ElectronScattering, "implementation of scattering by electrons")
        ATTRIBUTE_DEFAULT_VALUE(electronScattering, "Good")
        ATTRIBUTE_DISPLAYED_IF(electronScattering, "Level3")

        PROPERTY_BOOL(resonantScattering, "enable Lyman resonant scattering for all hydrogen-like ions")
        ATTRIBUTE_DEFAULT_VALUE(resonantScattering, "false")
        ATTRIBUTE_DISPLAYED_IF(resonantScattering, "Level2")
        ATTRIBUTE_RELEVANT_IF(includeThermalDispersion, "Lya")

    ITEM_END()

    //====================== Setup - Destruction =====================

public:
    ~XRayIonicGasMixFamily() override;

    void setupSelfBefore() override;

    //====================== Other functions ======================

public:
    vector<SnapshotParameter> parameterInfo() const override;

    const MaterialMix* mix(double Z, double T, const Array& parameters) override;

    const MaterialMix* mix() override;

private:
    // quick fix
    void setup();

    //======================== Data Members ========================

private:
    bool _setupDone{false};
    vector<string> _ionNames;
    XRayIonicGasMix::ElectronScattering _boundElectrons;
    vector<XRayIonicGasMix*> _mixes;
    XRayIonicGasMix* _defaultMix{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif

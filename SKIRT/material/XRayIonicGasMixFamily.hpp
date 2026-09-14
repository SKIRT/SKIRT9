
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIXFAMILY_HPP
#define XRAYIONICGASMIXFAMILY_HPP

#include "MaterialMixFamily.hpp"
#include "XRayIonicGasMix.hpp"

//////////////////////////////////////////////////////////////////////

/** The XRayIonicGasMixFamily class represents a family of XRayIonicGasMix instances, one for
    each spatial cell. This way the ion abundances and the temperature can vary across cells. The
    family defines a few properties that are shared by all mixes across all cells: \em ions, \em
    electronScattering, and \em resonantScattering. These properties can be found in the
	XRayIonicGasMix class.

    The \em ions property defines which ions must be imported. An ion can still have zero
    abundance in a given cell. As described in the XRayIonicGasMix class, that ion is then simply
    left out of the mix constructed for that cell.

    This family will also reuse mixes for cells with duplicate import parameters (ion abundances
	and temperature).
	
	This MaterialMixFamily is very memory-heavy 
	*/
class XRayIonicGasMixFamily : public MaterialMixFamily
{
    ENUM_DEF(ElectronScattering, None, Free, FreeWithPolarization, FreeBound)
        ENUM_VAL(ElectronScattering, None, "ignore electron")
        ENUM_VAL(ElectronScattering, Free, "use free-electron Compton scattering for all electrons")
        ENUM_VAL(ElectronScattering, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
        ENUM_VAL(ElectronScattering, FreeBound, "use an interpolation of free- and bound-electron Compton scattering")
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
    /** The destructor destructs all the XRayIonicGasMix instances that have been created. */
    ~XRayIonicGasMixFamily() override;

    /** This function calls the setup if not already done. This involves parsing the ions
        and electronScattering property, and adding a default mix. */
    void setupSelfBefore() override;

    //====================== Other functions ======================

public:
    /** This function returns the number and type of parameters used by family. For this class,
        this is all the relative abundances for the user-specified ions. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns (a pointer to) the XRayIonicGasMix corresponding to the given parameter
        values. The material mix family retains ownership of the returned material mix, and
        guarantees that it will not be destroyed until the family itself is destroyed.

        The number and type of parameters must match the information returned by the
        parameterInfo() function; if not the behavior is undefined. */
    const MaterialMix* mix(double Z, double T, const Array& parameters) override;

    /** This function returns (a pointer to) the default XRayIonicGasMix, corresponding to
        the parameters all set to zero. */
    const MaterialMix* mix() override;

private:
    /** This function performs the setup by parsing the user-configured ions and
        electronScattering property, and adding a default mix.*/
    void setup();

    //======================== Data Members ========================

private:
    bool _setupDone{false};
    vector<string> _ionNames;                             // parsed ion names
    XRayIonicGasMix::ElectronScattering _boundElectrons;  // parsed electronScattering property
    vector<XRayIonicGasMix*> _mixes;                      // all stored mixes
};

////////////////////////////////////////////////////////////////////

#endif

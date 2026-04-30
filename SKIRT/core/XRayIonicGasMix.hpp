/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XRAYIONICGASMIX_HPP
#define XRAYIONICGASMIX_HPP

#include "ArrayTable.hpp"
#include "DipolePhaseFunction.hpp"
#include "MaterialMix.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

class XRayIonicGasMix : public MaterialMix
{
    ENUM_DEF(ElectronScattering, None, Free, FreeWithPolarization, Good, Exact)
        ENUM_VAL(ElectronScattering, None, "ignore electron")
        ENUM_VAL(ElectronScattering, Free, "use free-electron Compton scattering for all electrons")
        ENUM_VAL(ElectronScattering, FreeWithPolarization,
                 "use free-electron Compton scattering with support for polarization")
        ENUM_VAL(ElectronScattering, Good, "use smooth Rayleigh scattering and exact bound-Compton scattering")
        ENUM_VAL(ElectronScattering, Exact, "use anomalous Rayleigh scattering and exact bound-Compton scattering")
    ENUM_END()

    ITEM_CONCRETE(XRayIonicGasMix, MaterialMix, "Ionised gas mix")
        ATTRIBUTE_TYPE_INSERT(XRayIonicGasMix, "GasMix,CustomMediumState")

        PROPERTY_STRING(ions, "the names of the ions for each element (e.g. H,He+,Li+1,..)")

        PROPERTY_DOUBLE_LIST(abundances, "the abundances of the ions in the same order as the ions property")

        PROPERTY_DOUBLE(temperature, "the temperature of the gas in K")
        ATTRIBUTE_QUANTITY(temperature, "temperature")
        ATTRIBUTE_MIN_VALUE(temperature, "[3")
        ATTRIBUTE_MAX_VALUE(temperature, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(temperature, "1e4")
        ATTRIBUTE_DISPLAYED_IF(temperature, "Level2")

        PROPERTY_ENUM(electronScattering, ElectronScattering, "implementation of scattering by electrons")
        ATTRIBUTE_DEFAULT_VALUE(electronScattering, "Free")
        ATTRIBUTE_DISPLAYED_IF(electronScattering, "Level3")

        PROPERTY_BOOL(resonantScattering, "enable Lyman resonant scattering for all hydrogen-like ions")
        ATTRIBUTE_DEFAULT_VALUE(resonantScattering, "false")
        ATTRIBUTE_DISPLAYED_IF(resonantScattering, "Level2")
        ATTRIBUTE_RELEVANT_IF(includeThermalDispersion, "Lya")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    explicit XRayIonicGasMix(SimulationItem* parent, string ions, vector<double> abundances, double temperature,
                             ElectronScattering boundElectrons, bool resonantScattering, bool setup);

    void setupSelfBefore() override;

    ~XRayIonicGasMix();

    //======== Private support functions =======

private:
    /** This function returns the index in the private wavelength grid corresponding to the
    specified wavelength. The parameters for converting a wavelength to the appropriate index
    are stored in data members during setup. */
    int indexForLambda(double lambda) const;

    double vtherm(int Z) const;

    //============= Capabilities =============

public:
    MaterialType materialType() const override;

    bool hasPolarizedScattering() const override;

    bool hasResonantScattering() const override;

    bool hasScatteringDispersion() const override;

    bool scatteringEmulatesSecondaryEmission() const override;

    bool hasLineEmission() const override;

    //============= Medium state setup =============

    vector<StateVariable> specificStateVariableInfo() const override;

    //============= Low-level material properties =============

    double mass() const override;

    double sectionAbs(double lambda) const override;

    double sectionSca(double lambda) const override;

    double sectionExt(double lambda) const override;

    //============= High-level photon life cycle =============

    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    Array sigmaSca(double lambda, const MaterialState* state) const;

    void setScatteringInfoIfNeeded(PhotonPacket* pp, const MaterialState* state, const double lambda) const;

    bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Temperature =======

public:
    /** This function returns an indicative temperature of the material mix when it would be
    embedded in a given radiation field. The implementation in this class ignores the radiation
    field and returns the (spatially constant) temperature configured for this material mix. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================

public:
    // base class for bound1-electron scattering helpers (public because we derive from it in anonymous namespace)
    class ScatteringHelper;

private:
    // all data is calculated in the setupSelfBefore(), but is persistent for use after setup

    struct IonParam
    {
        IonParam(short Z, short N) : Z(Z), N(N) {}

        short Z;  // atomic number
        short N;  // number of electrons
    };
    // Rayleigh scattering -> IonParam
    // Compton scattering -> IonParam
    // Fluorescence (+Lyman RC) -> FluorescenceParam
    struct FluorescenceParam
    {
        unsigned char Z;  // atomic number
        double lambda;    // wavelength (m)
        double width;     // width (eV)
    };
    // Resonant scattering -> LymanParam
    struct LymanParam
    {
        unsigned char Z;      // atomic number
        unsigned char index;  // Lyman index (alpha1/2, alpha3/2, beta1/2, ...)
        double lambda;        // wavelength (m)
        double a;             // Voigt parameter
        Array cumbranchingv;  // normalized cumulative branching
    };

    int _numIons;  // number of ions
    int _numFluo;  // number of fluorescence (+Lyman RC) transitions
    int _numLym;   // number of Lyman resonant scattering transitions

    // persistent data for scattering
    vector<IonParam> _ionParamv;
    vector<FluorescenceParam> _fluorescenceParamv;
    vector<LymanParam> _lymanParamv;
    vector<double> _vthermv;  // indexed on Z

    // wavelength grid (shifted to the left of the actually sampled points to approximate rounding)
    Array _lambdav;
    // cross sections
    Array _sigmaextv;              // indexed on lambdav
    Array _sigmascav;              // indexed on lambdav
    ArrayTable<2> _cumsigmascavv;  // indexed on lambdav, interactions (2*numIons + numFluo + numRes)

    // bound-electron scattering helpers depending on the configured implementation
    ScatteringHelper* _ray{nullptr};  // Rayleigh scattering helper
    ScatteringHelper* _com{nullptr};  // Compton scattering helper
    // dipole phase function for resonant scattering
    DipolePhaseFunction* _dpf{nullptr};
};

#endif

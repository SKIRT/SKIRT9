/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TRIVIALGASMIX_HPP
#define TRIVIALGASMIX_HPP

#include "MaterialMix.hpp"

////////////////////////////////////////////////////////////////////

/** The TrivialGasMix class represents a trivial gas mix defined by optical properties fully
    specified inside the configuration file. It is provided for testing purposes, and specifically
    for testing media with a negative absorption cross section (caused by stimulated emission).

    The class is designed for use in monochromatic simulations. The gas mix properties do not
    depend on wavelength, nor on any other property of the incoming photon packet or on any medium
    state variable (other than the number density). The gas mix supports absorption and
    Henyey-Greenstein scattering but no secondary emission. Each interacting entity in the gas is
    considered to be a hydrogen atom, which sets the conversion scale from mass to number
    normalization. If the gas distribution is normalized to a given optical depth, the values of
    the cross sections (see below) are essentially scale free.

    The following spatially-constant and wavelength-independent properties can be configured: the
    absorption cross section per hydrogen atom, which can be negative; the scattering cross section
    per hydrogen atom, which must be positive; and the scattering asymmetry parameter used with the
    Henyey-Greenstein phase function. */
class TrivialGasMix : public MaterialMix
{
    ITEM_CONCRETE(TrivialGasMix, MaterialMix, "A trivial gas mix for testing purposes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TrivialGasMix, "Level3")
        ATTRIBUTE_TYPE_INSERT(TrivialGasMix, "GasMix")

        PROPERTY_DOUBLE(absorptionCrossSection, "the absorption cross section per hydrogen atom")
        ATTRIBUTE_QUANTITY(absorptionCrossSection, "section")

        PROPERTY_DOUBLE(scatteringCrossSection, "the scattering cross section per hydrogen atom")
        ATTRIBUTE_QUANTITY(scatteringCrossSection, "section")
        ATTRIBUTE_MIN_VALUE(scatteringCrossSection, "[0")

        PROPERTY_DOUBLE(asymmetryParameter, "the scattering asymmetry parameter")
        ATTRIBUTE_MIN_VALUE(asymmetryParameter, "[-0.95")
        ATTRIBUTE_MAX_VALUE(asymmetryParameter, "0.95]")

    ITEM_END()

    //======== Capabilities =======

public:
    /** This function returns the fundamental material type represented by this material mix, which
        is MaterialType::Gas. */
    MaterialType materialType() const override;

    /** This function returns true if the configured extinction cross section (the sum of the
        configured absorption and scattering cross section) for this material mix is negative, and
        false otherwise. */
    bool hasNegativeExtinction() const override;

    //======== Medium state setup =======

public:
    /** This function returns a list of StateVariable objects describing the specific state
        variables used by the receiving material mix. For this class, the function returns just the
        descriptor for the number density. */
    vector<StateVariable> specificStateVariableInfo() const override;

    //======== Low-level material properties =======

public:
    /** This function returns the mass of a hydrogen atom. */
    double mass() const override;

    /** This function returns the absorption cross section per hydrogen atom configured for this
        material mix. The wavelength is not used. */
    double sectionAbs(double lambda) const override;

    /** This function returns the scattering cross section per hydrogen atom configured for this
        material mix. The wavelength is not used. */
    double sectionSca(double lambda) const override;

    /** This function returns the extinction cross section per hydrogen atom, i.e. the sum of the
        absorption and scattering cross sections per hydrogen atom configured for this material
        mix. The wavelength is not used. */
    double sectionExt(double lambda) const override;

    /** This function returns the scattering asymmetry parameter configured for this material mix.
        The wavelength is not used. */
    double asymmpar(double lambda) const override;

    //======== High-level photon life cycle =======

public:
    /** This function returns the absorption opacity for the number density given in the specified
        material state. The wavelength and the photon packet properties are not used. */
    double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the scattering opacity for the number density given in the specified
        material state. The wavelength and the photon packet properties are not used. */
    double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function returns the extinction opacity for the number density given in the specified
        material state. The wavelength and the photon packet properties are not used. */
    double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity for the given observer direction. The
        material state, wavelength and photon packet properties other than direction are not used.
        */
    void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs, Direction bfky,
                           const MaterialState* state, const PhotonPacket* pp) const override;

    /** This function performs a scattering event for this material mix on the specified photon
        packet. The photon packet properties other than luminosity and direction remain unchanged.
        */
    void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const override;

    //======== Temperature =======

public:
    /** This function returns an indicative temperature of zero. The specified material state and
        radiation field are noto used. */
    double indicativeTemperature(const MaterialState* state, const Array& Jv) const override;

    //======================== Data Members ========================
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SINGLEGRAINDUSTMIX_HPP
#define SINGLEGRAINDUSTMIX_HPP

#include "MaterialMix.hpp"
#include "DustExtinctionMixInterface.hpp"
#include "ScatteringMixInterface.hpp"
#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

/** SingleGrainDustMix is an abstract class describing a dust mix described by a single
    representative grain, with or without support for polarization by scattering. This base class
    includes the implementations of the required functions. Subclasses must merely provide the
    names of the relevant resource files. */
class SingleGrainDustMix : public MaterialMix, public DustExtinctionMixInterface, public ScatteringMixInterface
{
    ITEM_ABSTRACT(SingleGrainDustMix, MaterialMix, "a dust mix described by a single representative grain")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function opens the stored table resource(s) tabulating the representative grain
        properties. It invokes the resourceNameXxx() functions provided by each subclass to obtain
        the appropriate resource names. */
    void setupSelfBefore() override;

    /** This function must be implemented in each subclass. It returns the name of the stored table
        resource tabulating the basic optical properties (cross sections and asymmetry parameter)
        as a function of wavelength. */
    virtual string resourceNameForOpticalProps() const = 0;

    /** This function must be implemented in each subclass. To enable support for polarization by
        scattering, the function returns the name of the stored table resource tabulating the
        elements of the Mueller matrix as a function of wavelength and scattering angle. To disable
        polarization support, the function returns the empty string. */
    virtual string resourceNameForMuellerMatrix() const = 0;

    //============= Implementing DustExtinctionMixInterface =============

public:
    /** This function returns the absorption cross section per hydrogen atom
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sigmaabs(double lambda) const override;

    /** This function returns the scattering cross section per hydrogen atom
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sigmasca(double lambda) const override;

    /** This function returns the total extinction cross section per hydrogen atom
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ of the dust mix at wavelength \f$\lambda\f$. */
    double sigmaext(double lambda) const override;

    /** This function returns the absorption coefficient \f$\kappa^{\text{abs}}_\lambda\f$ of the
        dust mix at wavelength \f$\lambda\f$. */
    double kappaabs(double lambda) const override;

    /** This function returns the scattering coefficient \f$\kappa^{\text{sca}}_\lambda\f$ of the
        dust mix at wavelength \f$\lambda\f$. */
    double kappasca(double lambda) const override;

    /** This function returns the total extinction coefficient \f$\kappa^{\text{ext}}_\lambda =
        \kappa^{\text{abs}}_\lambda + \kappa^{\text{sca}}_\lambda\f$ of the dust mix at wavelength
        \f$\lambda\f$. */
    double kappaext(double lambda) const override;

    /** This function returns the scattering albedo \f$\varpi_\lambda =
        \varsigma_{\lambda}^{\text{sca}} / \varsigma_{\lambda}^{\text{ext}} =
        \kappa_{\lambda}^{\text{sca}} / \kappa_{\lambda}^{\text{ext}}\f$ of the dust mix at
        wavelength \f$\lambda\f$. */
    double albedo(double lambda) const override;

    //============= Implementing ScatteringMixInterface =============

public:
    /** This function returns true if this material mix supports polarization by scattering; false
        otherwise. */
    bool hasScatteringPolarization() const override;

    /** This function generates a new direction \f${\bf{k}}_{\text{new}}\f$ in case the specified
        photon packet scatters, and calculates the new polarization state of the scattered photon
        packet. The function passes the new direction to the caller as its return value, and stores
        the new polarization state in the provided Stokes vector. It is permitted for the provided
        Stokes vector to actually reside in the specified photon packet.

        For a dust mix that doesn't support polarization, the function generates the new direction
        from the normalized two-dimensional probability distribution \f[ p({\bf{k}}_{\text{new}})\,
        {\text{d}}{\bf{k}}_{\text{new}} = \Phi_\lambda({\bf{k}}_{\text{new}},
        {\bf{k}}_{\text{pp}})\, {\text{d}}{\bf{k}}_{\text{new}} \f] at the wavelength \f$\lambda\f$
        of the photon packet. Also, in that case, the provided Stokes vector is not modified.

        For a dust mix that does support polarization, the function generates a new direction
        \f${\bf{k}}_{\text{new}}\f$ after a scattering event, given that the original direction
        before the scattering event is \f${\bf{k}}\f$ and taking into account the polarization
        state of the photon. First, the polarization degree and angle are computed from the Stokes
        parameters. Then, scattering angles \f$\theta\f$ and \f$\phi\f$ are sampled from the phase
        function, and the Stokes vector is rotated into the scattering plane and transformed by
        applying the Mueller matrix. Finally, the new direction is computed from the previously
        sampled \f$\theta\f$ and \f$\phi\f$ angles. */
    Direction scatteringDirectionAndPolarization(StokesVector* out, const PhotonPacket* pp) const override;

    /** This function calculates the polarization state appropriate for a peel off photon packet
        generated by a scattering event for the specified photon packet, and stores the result in
        the provided Stokes vector.

        For a dust mix that doesn't support polarization, the provided Stokes vector is not
        modified, i.e. the function does nothing.

        For a dustmix that does support polarization, the function rotates the Stokes vector from
        the reference direction in the previous scattering plane into the peel-off scattering
        plane, applies the Mueller matrix on the Stokes vector, and further rotates the Stokes
        vector from the reference direction in the peel-off scattering plane to the x-axis of the
        instrument to which the peel-off photon packet is headed. */
    void scatteringPeelOffPolarization(StokesVector* out, const PhotonPacket* pp, Direction bfknew,
                                       Direction bfkx, Direction bfky) override;

    /** This function returns the value of the scattering phase function in case the specified
        photon packet is scattered to the specified new direction, where the phase function is
        normalized as \f[\int\Phi_\lambda(\Omega)\,\mathrm{d}\Omega=4\pi.\f]

        For a dustmix that doesn't support polarization, the function returns
        \f$\Phi_\lambda({\bf{k}}_{\text{pp}}, {\bf{k}}_{\text{new}})\f$ for the current propagation
        direction of the photon packet \f${\bf{k}}_{\text{pp}}\f$ and the specified new direction
        \f${\bf{k}}_{\text{new}}\f$, at the wavelength \f$\lambda\f$ of the photon packet.
        Specifically, the phase function is taken to be the Henyey-Greenstein phase function
        parameterized by the asymmetry parameter \f$g_\lambda\f$ of the dust mix at the wavelength
        \f$\lambda\f$ of the photon packet.

        For a dustmix that does support polarization, the function returns the phase function for
        polarized radiation given by \f[\Phi_\lambda(\Omega) = N \left( S_{11,\lambda}(\theta) +
        P_\text{L} S_{12,\lambda}(\theta) \cos 2(\varphi-\gamma) \right)\f] where \f$\lambda\f$ is
        the wavelength of the photon packet, \f$\theta\f$ is the angle between the photon packet's
        propagation direction and the new scattering direction; \f$\phi\f$ is the angle between the
        previous and current scattering plane of the photon packet; \f$\gamma\f$ is the
        polarization angle of the photon packet, \f$P_\text{L}\f$ is the linear polarization degree
        of the photon packet; and \f$N\f$ is a normalization factor to ensure that the integral
        over the unit sphere is equal to \f$4\pi\f$. */
    double phaseFunctionValue(const PhotonPacket* pp, Direction bfknew) const override;

private:
    /** This function returns a random scattering angle \f$\theta\f$ sampled from the phase
        function for a given wavelength \f$\lambda\f$. */
    double sampleTheta(double lambda) const;

    /** This function returns a random scattering angle \f$\phi\f$ sampled from the phase function
        according to the scattering angle \f$\theta\f$ and the incident linear polarization degree
        and polarization angle, at wavelength \f$\lambda\f$. */
    double samplePhi(double lambda, double theta, double polDegree, double polAngle) const;

    //======================== Data Members ========================

private:
    // basic properties - initialized during setup
    StoredTable<1> _sigmaabs;       // indexed on lambda
    StoredTable<1> _sigmasca;       // indexed on lambda
    StoredTable<1> _asymmpar;       // indexed on lambda
    double _mu{0};

    // polarization properties - initialized during setup
    bool _polarization{false};
    StoredTable<2> _S11;            // indexed on lambda, theta
    StoredTable<2> _S12;            // indexed on lambda, theta
    StoredTable<2> _S33;            // indexed on lambda, theta
    StoredTable<2> _S34;            // indexed on lambda, theta

    // precalculated discretizations - initialized during setup
    Array _thetav;                  // indexed on t
    Array _phiv;                    // indexed on f
    Array _phi1v;                   // indexed on f
    Array _phisv;                   // indexed on f
    Array _phicv;                   // indexed on f
};

////////////////////////////////////////////////////////////////////

#endif

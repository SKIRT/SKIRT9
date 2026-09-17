/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ELECTRONSCATTERINGHELPER_HPP
#define ELECTRONSCATTERINGHELPER_HPP

#include "ComptonPhaseFunction.hpp"
#include "DipolePhaseFunction.hpp"
#include "Direction.hpp"
#include "Range.hpp"
#include "SimulationItem.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

/** This namespace offers a set of helper classes that implement the various treatments of photon
    scattering by electrons -- bound to an atom or ion, or free -- supported by the XRayAtomicGasMix
    and XRayIonicGasMix material mix classes. Both classes select one concrete helper for Rayleigh
    (elastic) scattering and one for Compton (inelastic) scattering, based on their own configured
    implementation options, and delegate all scattering-related calculations to these two helpers
    for the remainder of the simulation; see the documentation of those classes for the physics
    (cross section and phase function formulas, energy shift) implemented by each helper.

    All helpers derive from the abstract base class Helper, which declares the common interface:
    a scattering cross section per atom or ion, and peel-off and perform-scattering functions with
    and without support for polarization. A helper that does not (or does) support polarization
    needs to override only the unpolarized (or only the polarized) pair of functions; Helper
    provides a default implementation for the other pair (a no-op for the unpolarized functions,
    and a fall-through to the unpolarized functions for the polarized ones) so that a derived class
    never needs to implement both. All functions other than the cross section take both the atomic
    number \f$Z\f$ and the number of bound electrons \f$N\f$ of the scattering species; because
    XRayAtomicGasMix always represents neutral atoms (\f$N=Z\f$), Helper also offers a convenience
    overload of each function that takes only \f$Z\f$.

    The concrete helpers are:

    - NoScatteringHelper: a stub with zero cross section, used when a scattering channel (Rayleigh
    or Compton) is disabled altogether.

    - FreeComptonHelper and FreeComptonWithPolarizationHelper: Compton scattering by free electrons,
    delegating the cross section, energy shift and phase function to the ComptonPhaseFunction class
    (using the Thomson limit at low energies, where the Compton formulas become numerically
    unstable), without or with support for polarization, respectively.

    - BoundComptonHelper: Compton scattering by electrons bound to an atom, using tabulated cross
    sections, incoherent scattering functions and target-electron momentum distributions.

    - FreeBoundComptonHelper: a weighted combination of free- and bound-electron Compton scattering
    for an ion with \f$N\f$ out of \f$Z\f$ electrons still bound, used only by XRayIonicGasMix.

    - SmoothRayleighHelper and AnomalousRayleighHelper: two levels of fidelity for Rayleigh
    (coherent, elastic) scattering, using tabulated cross sections and atomic form factors, with
    AnomalousRayleighHelper additionally including tabulated anomalous (energy-dependent) scattering
    corrections; both fall back to \f$Z^2\f$ Thomson scattering below the energy range of the
    tabulated data.

    Aside from NoScatteringHelper and FreeBoundComptonHelper, all of these helpers load their
    tabulated data, indexed on atomic number up to \f$Z=30\f$, from resource files during
    construction. */
namespace ElectronScatteringHelper
{
    /** This is the abstract base class for the electron-scattering helpers offered by this
        namespace; see the namespace documentation for the overall design and the list of concrete
        helpers. */
    class Helper
    {
    public:
        virtual ~Helper();

        // return scattering cross section for atom in m2
        virtual double sectionSca(double lambda, int Z, int N) const = 0;

        // convenience overload for a neutral atom, i.e. N=Z
        double sectionSca(double lambda, int Z) const { return sectionSca(lambda, Z, Z); }

        // peel-off unpolarized scattering event: override this in helpers that don't support polarization
        virtual void peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk, Direction bfkobs) const;

        // convenience overload for a neutral atom, i.e. N=Z
        void peeloffScattering(double& I, double& lambda, int Z, Direction bfk, Direction bfkobs) const
        {
            peeloffScattering(I, lambda, Z, Z, bfk, bfkobs);
        }

        // perform unpolarized scattering event: override this in helpers that don't support polarization
        virtual Direction performScattering(double& lambda, int Z, int N, Direction bfk) const;

        // convenience overload for a neutral atom, i.e. N=Z
        Direction performScattering(double& lambda, int Z, Direction bfk) const
        {
            return performScattering(lambda, Z, Z, bfk);
        }

        // peel-off polarized scattering event: override this in helpers that do support polarization
        virtual void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, int Z, int N,
                                       Direction bfk, Direction bfkobs, Direction bfky, const StokesVector* sv) const;

        // convenience overload for a neutral atom, i.e. N=Z
        void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, int Z, Direction bfk,
                               Direction bfkobs, Direction bfky, const StokesVector* sv) const
        {
            peeloffScattering(I, Q, U, V, lambda, Z, Z, bfk, bfkobs, bfky, sv);
        }

        // perform polarized scattering event: override this in helpers that do support polarization
        virtual Direction performScattering(double& lambda, int Z, int N, Direction bfk, StokesVector* sv) const;

        // convenience overload for a neutral atom, i.e. N=Z
        Direction performScattering(double& lambda, int Z, Direction bfk, StokesVector* sv) const
        {
            return performScattering(lambda, Z, Z, bfk, sv);
        }
    };

    ////////////////////////////////////////////////////////////////////

    /** This class implements a stub scattering helper with a zero cross section, used when a given
        type of scattering (Rayleigh or Compton) is disabled altogether. */
    class NoScatteringHelper : public Helper
    {
    public:
        NoScatteringHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;
    };

    ////////////////////////////////////////////////////////////////////

    /** This class implements Compton scattering by free electrons, without support for
        polarization. The cross section, phase function and energy shift are delegated to the
        ComptonPhaseFunction class, except at low energies (below 0.1 keV), where the Compton
        formulas become numerically unstable and the Thomson limit is used instead. */
    class FreeComptonHelper : public Helper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;

        void peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk, Direction bfkobs) const override;

        Direction performScattering(double& lambda, int Z, int N, Direction bfk) const override;
    };

    ////////////////////////////////////////////////////////////////////

    /** This class implements Compton scattering by free electrons, with support for polarization.
        As for FreeComptonHelper, the cross section, phase function and energy shift are delegated
        to the ComptonPhaseFunction class, except at low energies (below 0.1 keV), where the
        Thomson limit is used instead. */
    class FreeComptonWithPolarizationHelper : public Helper
    {
    private:
        ComptonPhaseFunction _cpf;
        DipolePhaseFunction _dpf;

    public:
        FreeComptonWithPolarizationHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;

        void peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, int Z, int N, Direction bfk,
                               Direction bfkobs, Direction bfky, const StokesVector* sv) const override;

        Direction performScattering(double& lambda, int Z, int N, Direction bfk, StokesVector* sv) const override;
    };

    ////////////////////////////////////////////////////////////////////

    /** This class implements Compton scattering by electrons bound to an atom, without support for
        polarization, using tabulated cross sections, incoherent scattering functions, and
        target-electron momentum distributions (indexed on atomic number) to calculate the cross
        section, phase function, and the resulting photon energy shift. */
    class BoundComptonHelper : public Helper
    {
    private:
        // resources loaded from file
        vector<Array> _CSv;  // 0: E (keV->1); 1-30: bound Compton cross sections (cm2->m2)
        vector<Array> _SFv;  // 0: q (1); 1-30: incoherent scattering functions (1)
        vector<Array> _CPv;  // 0: E (keV->1); 1-30: pdf for target electron momentum (1)
        vector<Array> _IBv;  // 0: E (keV->1) ionisation energy of the outer subshell electrons

        // precalculated cumulative distributions for target electron momentum
        Range _cumRange;
        vector<Array> _cumCPv;  // 0: E axis; 1-30: cumulative pdf for target electron momentum

        // precalculated discretizations
        Array _costhetav;
        Array _sinthetav;
        Array _sin2thetav;
        Array _sintheta2v;

        // cache
        Random* _random{nullptr};

    public:
        BoundComptonHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const;

        double generateCosineFromPhaseFunction(double x, double Z) const;

        // sample a target electron momentum from the distribution with the given maximum
        double sampleMomentum(double pmax, double Z) const;

        // returns the augmented inverse Compton factor
        double augmentedInverseComptonFactor(double x, double costheta, double Z) const;

    public:
        void peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk, Direction bfkobs) const override;

        Direction performScattering(double& lambda, int Z, int N, Direction bfk) const override;
    };

    ////////////////////////////////////////////////////////////////////

    /** This class implements Compton scattering for an ion with \f$N\f$ out of \f$Z\f$ electrons
        still bound, without support for polarization, as a weighted combination of free-electron
        (FreeComptonHelper) and bound-electron (BoundComptonHelper) Compton scattering, with weights
        \f$1-N/Z\f$ and \f$N/Z\f$, respectively. This helper is used only by XRayIonicGasMix, which
        represents ions rather than neutral atoms. */
    class FreeBoundComptonHelper : public Helper
    {
    private:
        FreeComptonHelper _free;
        BoundComptonHelper _bound;

        Random* _random{nullptr};

    public:
        FreeBoundComptonHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;

    public:
        void peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk, Direction bfkobs) const override;

        Direction performScattering(double& lambda, int Z, int N, Direction bfk) const override;
    };

    ////////////////////////////////////////////////////////////////////

    /** This class implements Rayleigh (coherent, elastic) scattering, without support for
        polarization, using tabulated cross sections and atomic form factors (indexed on atomic
        number) to calculate the cross section and phase function. Below the energy range of the
        tabulated data, \f$Z^2\f$ Thomson scattering is used instead. */
    class SmoothRayleighHelper : public Helper
    {
    private:
        vector<Array> _RSSv;  // 0: E (keV->1); 1-30: smooth Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-30: atomic form factors (1)
        Random* _random{nullptr};
        DipolePhaseFunction _dpf;

        // precalculated discretizations
        Array _costhetav;
        Array _cos2thetav;
        Array _sinthetav;
        Array _sintheta2v;

    public:
        SmoothRayleighHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const;

        double generateCosineFromPhaseFunction(double x, double Z) const;

    public:
        void peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk, Direction bfkobs) const override;

        Direction performScattering(double& lambda, int Z, int N, Direction bfk) const override;
    };

    ////////////////////////////////////////////////////////////////////

    /** This class implements Rayleigh (coherent, elastic) scattering, without support for
        polarization, at a higher level of fidelity than SmoothRayleighHelper: in addition to
        tabulated cross sections and atomic form factors, it uses tabulated real and imaginary
        anomalous scattering functions (indexed on atomic number) that account for the
        energy-dependent deviation from the smooth form-factor approximation near absorption edges.
        Below the energy range of the tabulated data, \f$Z^2\f$ Thomson scattering is used instead.
        */
    class AnomalousRayleighHelper : public Helper
    {
    private:
        vector<Array> _RSAv;  // 2*Z: E (keV->1); 2*Z+1: anomalous Rayleigh cross sections (cm2->m2)
        vector<Array> _FFv;   // 0: q (1); 1-30: atomic form factors (1)
        vector<Array> _F1v;   // 2*Z: E (keV->1); 2*Z+1: Real anomalous scattering function (1)
        vector<Array> _F2v;   // 2*Z: E (keV->1); 2*Z+1: Imaginary anomalous scattering function (1)
        Random* _random{nullptr};
        DipolePhaseFunction _dpf;

        // precalculated discretizations
        Array _costhetav;
        Array _cos2thetav;
        Array _sinthetav;
        Array _sintheta2v;

    public:
        AnomalousRayleighHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;

    private:
        double phaseFunctionValue(double x, double costheta, int Z) const;

        double generateCosineFromPhaseFunction(double x, double Z) const;

    public:
        void peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk, Direction bfkobs) const override;

        Direction performScattering(double& lambda, int Z, int N, Direction bfk) const override;
    };
}

////////////////////////////////////////////////////////////////////

#endif

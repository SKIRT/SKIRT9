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

namespace ElectronScatteringHelper
{
    // ---- base class for scattering helpers ----

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

    // ---- no scattering helper ----

    // this helper does nothing; it is used as a stub in case there is no scattering of a given type
    class NoScatteringHelper : public Helper
    {
    public:
        NoScatteringHelper(SimulationItem* item);

        double sectionSca(double lambda, int Z, int N) const override;
    };

    ////////////////////////////////////////////////////////////////////

    // ---- free-electron Compton scattering helper ----

    // this helper forwards all calls to an external helper class for regular Compton scattering
    // (or Thomson scattering for lower energies, because Compton becomes numerically unstable)
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

    // ---- free-electron Compton with polarization scattering helper ----

    // this helper forwards all calls to an external helper class for Compton scattering
    // (or Thomson scattering for lower energies) with support for polarization
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

    // ---- bound-electron Compton scattering helper ----

    // this helper implements bound-electron Compton scattering
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

    // ---- free-bound Compton scattering helper ----

    // this helper implements an interpolation of free- and bound-Compton scattering
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

    // ---- smooth Rayleigh scattering helper ----

    // this helper implements smooth Rayleigh scattering;
    // below the energy limit of the tabulated data, use Thomson scattering instead
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

    // ---- anomalous Rayleigh scattering helper ----

    // this helper implements anomalous Rayleigh scattering
    // below the energy limit of the tabulated data, use Thomson scattering instead
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
        double phaseFunctionValue(double x, double costheta, int Z, int N) const;

        double generateCosineFromPhaseFunction(double x, double Z) const;

    public:
        void peeloffScattering(double& I, double& lambda, int Z, int N, Direction bfk, Direction bfkobs) const override;

        Direction performScattering(double& lambda, int Z, int N, Direction bfk) const override;
    };
}

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTOIONIZATIONRATES_HPP
#define PHOTOIONIZATIONRATES_HPP

#include "Basics.hpp"

//////////////////////////////////////////////////////////////////////

/** The PhotoIonizationRates class provides rate coefficients for recombination (radiative and
    dielectronic), collisional ionization, and charge exchange for all ions of H, He, C, N, O, Ne,
    Mg, Si, S, and Fe. Sources: Case B recombination from Hui & Gnedin (1997), radiative
    recombination from Badnell (2006), dielectronic recombination from Badnell (2006), collisional
    ionization from Voronov (1997), charge exchange from Kingdon & Ferland (1996).
    All rates in cm^3/s, temperature T in K.

    The individual per-ion rate functions and their fitting coefficients are file-local to the
    implementation translation unit; only the dispatch entry points and the small set of
    rates consumed directly by the solver / cooling code are exposed here. */
class PhotoIonizationRates final
{
public:
    // ============================================================
    //  Dispatch functions by linear ion index
    //  (same ordering as VernerCrossSections)
    // ============================================================

    /** Total recombination rate alpha_RR + alpha_DR for the recombination from ion stage m+1 -> m.
        The ion index refers to the LOWER state (the state being recombined INTO).
        For H and He, returns Case B recombination. For metals, returns RR + DR. */
    static double alphaTotal(int ion, double T);

    /** Collisional ionization rate [cm^3/s] for ion stage m -> m+1.
        The ion index refers to the state being ionized FROM. */
    static double beta(int ion, double T);

    /** Hydrogen charge exchange recombination rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HI -> ion_lower + HII.
        The ion index refers to the HIGHER stage (the one being recombined). */
    static double xiHI(int ion, double T);

    /** Hydrogen charge exchange ionization rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HII -> ion_higher + HI.
        The ion index refers to the LOWER stage (the one being ionized). */
    static double chiHII(int ion, double T);

    /** Helium charge exchange recombination rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HeI -> ion_lower + HeII.
        The ion index refers to the HIGHER stage (the one being recombined). */
    static double xiHeI(int ion, double T);

    /** Helium charge exchange ionization rate [cm^3/s] for a given ion (Verner index).
        Reaction: ion + HeII -> ion_higher + HeI.
        The ion index refers to the LOWER stage (the one being ionized). */
    static double chiHeII(int ion, double T);

    // ============================================================
    //  H and He rates used directly by the solver and cooling code
    // ============================================================

    /** Case B recombination coefficient for HII [cm^3/s]. Hui & Gnedin (1997). */
    static double alphaBHII(double T);

    /** Case B recombination coefficient for HeII [cm^3/s]. Hui & Gnedin (1997). */
    static double alphaBHeII(double T);

    /** Case B recombination coefficient for HeIII [cm^3/s]. Hui & Gnedin (1997). */
    static double alphaBHeIII(double T);

    /** Dielectronic recombination rate for HeII [cm^3/s]. Aldrovandi & Pequignot (1973). */
    static double alphaDRHeII(double T);

    /** HI collisional ionization rate [cm^3/s]. */
    static double betaHI(double T);

    /** HeI collisional ionization rate [cm^3/s]. */
    static double betaHeI(double T);

    /** HeII collisional ionization rate [cm^3/s]. */
    static double betaHeII(double T);

    /** HI charge exchange with HeII [cm^3/s]. */
    static double xiHIHeII(double T);

    /** HI charge exchange with HeIII [cm^3/s]. */
    static double xiHIHeIII(double T);
};

//////////////////////////////////////////////////////////////////////

#endif

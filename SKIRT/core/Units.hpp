/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef UNITS_HPP
#define UNITS_HPP

#include "SimulationItem.hpp"
#include "SkirtUnitDef.hpp"

//////////////////////////////////////////////////////////////////////

/** Units is the abstract base class of the classes that represent the unit systems supported by
    SKIRT for input/output purposes (internally in SKIRT, everything is in SI units). The Units
    class serves several purposes.

    Firstly, the Units base class and its derived classes enable the SMILE units mechanism in SKIRT
    parameter files, as described in the documentation of the UnitDef class.

    Secondly, the Units class allows the user to configure the flux output style as an attribute in
    the SKIRT parameter file.

    Finally, the Units class offers functionality for use by other classes in the simulation item
    hierarchy. It provides functions for converting physical quantities from internal SI units to
    external output units depending on the unit system (determined by the name of the subclass
    being called) and, where applicable, on the selected flux output style. */
class Units : public SimulationItem
{
    /** The enumeration type indicating the output style for flux density and surface brightness.
        Neutral indicates \f$\lambda F_\lambda = \nu F_\nu\f$; Wavelength indicates
        \f$F_\lambda\f$; and Frequency indicates \f$F_\nu\f$. */
    ENUM_DEF(FluxOutputStyle, Neutral, Wavelength, Frequency)
        ENUM_VAL(FluxOutputStyle, Neutral, "neutral: λ F_λ = ν F_ν")
        ENUM_VAL(FluxOutputStyle, Wavelength, "per unit of wavelength: F_λ")
        ENUM_VAL(FluxOutputStyle, Frequency, "per unit of frequency: F_ν")
    ENUM_END()

    ITEM_ABSTRACT(Units, SimulationItem, "a units system")

        PROPERTY_ENUM(fluxOutputStyle, FluxOutputStyle, "the output style for flux density and surface brightness")
        ATTRIBUTE_DEFAULT_VALUE(fluxOutputStyle, "Frequency")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns true if the specified combination of physical quantity and unit or
        unit system is present in the unit definition, and false if not. */
    bool has(string qty, string unit) const;

    /** This function returns a string containing the name of the output unit adopted by the
        program for the specified physical quantity. The name of the physical quantity must be
        specified in all lowercase and without any spaces. The function throws a fatal error if the
        specified physical quantity is unknown. */
    string unit(string qty) const;

    /** This function converts a physical value from the specified units to internal program units.
        If the specified combination is not present in the unit definition, the function throws an
        exception. */
    double in(string qty, string unit, double value) const;

    /** This function converts a physical value from internal SI units to the output units adopted
        by the program. The name of the physical quantity must be specified in all lowercase and
        without any spaces. If the specified combination of physical quantity and unit is unknown,
        the function throws a fatal error. */
    double out(string qty, double value) const;

    /** This function returns a string containing the name of the unit of length adopted by the
        program for output. Apart from this unit of length, the program uses an independent unit of
        distance, because usually the distance to an astronomical object and the typical length
        scale within an astronomical object are quite different. */
    string ulength() const;

    /** This function converts the length \f$x\f$ from the internally used SI units (m) to the
        program's output units. */
    double olength(double x) const;

    /** This function returns a string containing the name of the unit of distance adopted by the
        program for output. Apart from this unit of distance, the program uses an independent unit
        of length, because usually the distance to an astronomical object and the typical length
        scale within an astronomical object are quite different. */
    string udistance() const;

    /** This function converts the distance \f$d\f$ from the internally used SI units (m) to the
        program's output units. */
    double odistance(double d) const;

    /** This function returns a string containing the name of the unit of wavelength adopted by the
        program for output. */
    string uwavelength() const;

    /** This function converts the wavelength \f$\lambda\f$ from the internally used SI units (m)
        to the program's output units. */
    double owavelength(double lambda) const;

    /** This function returns a string containing the name of the unit of dust grain size adopted
        by the program for output. */
    string ugrainsize() const;

    /** This function converts the dust grain size \f$a\f$ from the internally used SI units (m) to
        the program's output units. */
    double ograinsize(double a) const;

    /** This function returns a string containing the name of the unit of cross section adopted by
        the program for output. */
    string usection() const;

    /** This function converts the volume \f$C\f$ from the internally used SI units
        (\f${\text{m}}^2\f$) to the program's output units. */
    double osection(double C) const;

    /** This function returns a string containing the name of the unit of volume adopted by the
        program for output. */
    string uvolume() const;

    /** This function converts the volume \f$V\f$ from the internally used SI units
        (\f${\text{m}}^3\f$) to the program's output units. */
    double ovolume(double V) const;

    /** This function returns a string containing the name of the unit of velocity adopted by the
        program for output. */
    string uvelocity() const;

    /** This function converts the velocity \f$v\f$ from the internally used SI units
        (\f${\text{m}}\, {\text{s}}^{-1}\f$) to the program's output units. */
    double ovelocity(double v) const;

    /** This function returns a string containing the name of the unit of mass adopted by the
        program for output. */
    string umass() const;

    /** This function converts the mass \f$M\f$ from the internally used SI units (kg) to the
        program's output units. */
    double omass(double M) const;

    /** This function returns a string containing the name of the unit of bulk mass adopted by the
        program for output. */
    string ubulkmass() const;

    /** This function converts the bulk mass \f$\mu\f$ from the internally used SI units (kg) to
        the program's output units. */
    double obulkmass(double mu) const;

    /** This function returns a string containing the name of the unit of bulk mass density adopted
        by the program for output. */
    string ubulkmassdensity() const;

    /** This function converts the bulk mass density \f$\rho\f$ from the internally used SI units
        (\f${\text{kg}}\, {\text{m}}^{-3}\f$) to the program's output units. */
    double obulkmassdensity(double rho) const;

    /** This function returns a string containing the name of the unit of mass surface density adopted
        by the program for output. */
    string umasssurfacedensity() const;

    /** This function converts the mass surface density \f$\Sigma\f$ from the internally used SI units
        (\f${\text{kg}}\, {\text{m}}^{-2}\f$) to the program's output units. */
    double omasssurfacedensity(double Sigma) const;

    /** This function returns a string containing the name of the unit of mass volume density adopted by
        the program for output. */
    string umassvolumedensity() const;

    /** This function converts the mass volume density \f$\rho\f$ from the internally used SI units
        (\f${\text{kg}}\, {\text{m}}^{-3}\f$) to the program's output units. */
    double omassvolumedensity(double rho) const;

    /** This function returns a string containing the name of the unit of number surface density adopted
        by the program for output. */
    string unumbersurfacedensity() const;

    /** This function converts the number surface density \f$\mathcal{N}\f$ from the internally used SI units
        (\f${\text{m}}^{-2}\f$) to the program's output units. */
    double onumbersurfacedensity(double N) const;

    /** This function returns a string containing the name of the unit of number volume density adopted by
        the program for output. */
    string unumbervolumedensity() const;

    /** This function converts the number volume density \f$n\f$ from the internally used SI units
        (\f${\text{m}}^{-3}\f$) to the program's output units. */
    double onumbervolumedensity(double n) const;

    /** This function returns a string containing the name of the unit of mass coefficient adopted
        by the program for output. */
    string umasscoefficient() const;

    /** This function converts the mass coefficient \f$\kappa\f$ from the internally used SI units
        (\f${\text{m}}^{-2}\, {\text{kg}}^{-1}\f$) to the program's output units. */
    double omasscoefficient(double kappa) const;

    /** This function returns a string containing the name of the unit of energy adopted by the
        program for output. */
    string uenergy() const;

    /** This function converts the energy \f$E\f$ from the internally used SI units (J) to the
        program's output units. */
    double oenergy(double E) const;

    /** This function returns a string containing the name of the unit of magnetic field adopted by
        the program for output. */
    string umagneticfield() const;

    /** This function converts the magnetic field strength \f$B\f$ from the internally used SI
        units (T) to the program's output units. */
    double omagneticfield(double B) const;

    /** This function returns a string containing the name of the unit of bolometric luminosity
        adopted by the program for output. */
    string ubolluminosity() const;

    /** This function converts the bolometric luminosity \f$L\f$ from the internally used SI units
        (W) to the program's output units. */
    double obolluminosity(double L) const;

    /** This function returns a string describing the monochromatic luminosity output style adopted
        by the program. */
    string smonluminosity() const;

    /** This function returns a string containing the name of the unit of monochromatic luminosity
        adopted by the program for output, depending on the selected flux output style. */
    string umonluminosity() const;

    /** This function converts the monochromatic per-wavelength luminosity \f$L_\lambda\f$ from the
        internally used SI units (\f${\text{W}}\, {\text{m}}^{-1}\f$) to the program's flux output
        style (neutral, wavelength or frequency) and units. */
    double omonluminosityWavelength(double lambda, double Llambda) const;

    /** This function converts the monochromatic per-frequency luminosity \f$L_\nu\f$ from the
        internally used SI units (\f${\text{W}}\, {\text{Hz}}^{-1}\f$) to the program's flux output
        style (neutral, wavelength or frequency) and units. */
    double omonluminosityFrequency(double lambda, double Lnu) const;

    /** This function returns a string describing the flux density output style adopted by the
        program. */
    string sfluxdensity() const;

    /** This function returns a string containing the name of the unit of flux density adopted by
        the program for output, depending on the selected flux output style. */
    string ufluxdensity() const;

    /** This function converts the per-wavelength flux density \f$F_\lambda\f$ for wavelength
        \f$\lambda\f$ from the internally used SI units (\f${\text{W}}\, {\text{m}}^{-3}\f$) to the
        program's flux output style (neutral, wavelength or frequency) and units. */
    double ofluxdensityWavelength(double lambda, double Flambda) const;

    /** This function converts the per-frequency flux density \f$F_\nu\f$ for wavelength
        \f$\lambda\f$ from the internally used SI units (\f${\text{W}}\, {\text{m}}^{-2}\,
        {\text{Hz}}^{-1}\f$) to the program's flux output style (neutral, wavelength or frequency)
        and units. */
    double ofluxdensityFrequency(double lambda, double Fnu) const;

    /** This function returns a string describing the surface brightness output style adopted by
        the program. */
    string ssurfacebrightness() const;

    /** This function returns a string containing the name of the unit of surface brightness
        adopted by the program for output, depending on the selected flux output style. */
    string usurfacebrightness() const;

    /** This function converts the per-wavelength surface brightness \f$f_\lambda\f$ for wavelength
        \f$\lambda\f$ from the internally used SI units (\f${\text{W}}\, {\text{m}}^{-3}\,
        {\text{sr}}^{-1}\f$) to the program's flux output style (neutral, wavelength or frequency)
        and units. */
    double osurfacebrightnessWavelength(double lambda, double flambda) const;

    /** This function converts the per-frequency surface brightness \f$f_\nu\f$ for wavelength
        \f$\lambda\f$ from the internally used SI units (\f${\text{W}}\, {\text{m}}^{-2}\,
        {\text{Hz}}^{-1}\, {\text{sr}}^{-1}\f$) to the program's flux output style (neutral,
        wavelength or frequency) and units. */
    double osurfacebrightnessFrequency(double lambda, double fnu) const;

    /** This function returns a string describing the mean intensity output style adopted by
        the program. */
    string smeanintensity() const;

    /** This function returns a string containing the name of the unit of mean intensity
        adopted by the program for output, depending on the selected flux output style. */
    string umeanintensity() const;

    /** This function converts the per-wavelength mean intensity \f$J_\lambda\f$ for wavelength
        \f$\lambda\f$ from the internally used SI units (\f${\text{W}}\, {\text{m}}^{-3}\,
        {\text{sr}}^{-1}\f$) to the program's flux output style (neutral, wavelength or frequency)
        and units. */
    double omeanintensityWavelength(double lambda, double Jlambda) const;

    /** This function converts the per-frequency mean intensity \f$J_\nu\f$ for wavelength
        \f$\lambda\f$ from the internally used SI units (\f${\text{W}}\, {\text{m}}^{-2}\,
        {\text{Hz}}^{-1}\, {\text{sr}}^{-1}\f$) to the program's flux output style (neutral,
        wavelength or frequency) and units. */
    double omeanintensityFrequency(double lambda, double Jnu) const;

    /** This function returns a string containing the name of the unit of temperature adopted by
        the program for output. */
    string utemperature() const;

    /** This function converts the temperature \f$T\f$ from the internally used SI units (K) to the
        program's output units. */
    double otemperature(double T) const;

    /** This function returns a string containing the name of the unit of angular size adopted by
        the program for output. */
    string uangle() const;

    /** This function converts the angular size \f$\theta\f$ from the internally used SI units
        (rad) to the program's output units. */
    double oangle(double theta) const;

    /** This function returns a string containing the name of the unit of positioning angle adopted
        by the program for output. */
    string uposangle() const;

    /** This function converts the positioning angle \f$\theta\f$ from the internally used SI units
        (rad) to the program's output units. */
    double oposangle(double theta) const;

    /** This function returns a string containing the name of the unit of solid angle adopted by
        the program for output. */
    string usolidangle() const;

    /** This function converts the solid angle \f$\Omega\f$ from the internally used SI units (sr)
        to the program's output units. */
    double osolidangle(double Omega) const;

    /** This function returns a string containing the name of the unit of pressure adopted by the
        program for output. */
    string upressure() const;

    /** This function converts the pressure \f$p\f$ from the internally used SI units (Pa) to the
        program's output units. */
    double opressure(double p) const;

    //======================== Data Members ========================

private:
    SkirtUnitDef _unitDef;
};

//////////////////////////////////////////////////////////////////////

#endif

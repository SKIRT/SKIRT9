/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSEDFAMILY_HPP
#define TODDLERSSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** TODDLERS = Time evolution of Dust Diagnostics and Line Emission from Regions
    containing young Stars. 
    
    An instance of the ToddlersSEDFamily class represents the family of 
    Toddlers star-region template SEDs, parameterized on age, metallicity,
    star formation efficiency, and cloud number density. The 
    SED is scaled by the particle mass.
    The scaling assumes a cloud mass function obeying a power law with slope of -1.8, i.e.,
    dN/dM \propto M^{-1.8} running from cloud mass of 1e5 - 1e6.75 Msun.
    The model spectra have 9756 wavelengths and cover 90 ages from .1-30Myr, 5 metallicity values from
    .001 to .04 (solar being .014), 7 star formation efficiencies from 1 - 15 %, and 9 cloud number
    densities ranging from 10 cm^{-3} to 2560 cm^{-3}.
	The PAH fraction is a maximal value scaled by the neutral Hydrogen abundance in a shell, 
	two maximum q_PAH values are avaiable: High(Default) = 4.6 %, Low = 1 %.
    Additionally, this class allows the use of low resolution (R=300) as well as high resolution spectra
    to represent the total SED (line + continuum). In the case of high resolution case, the lines are 
    at high resolution while the continua is the same as the low resolution case.
    The SEDs are tabulated over a wavelength range from 100 Angstrom to 3000 micron.
    The data was generated using the Cloudy master branch.
    
    For more details refer to Kapoor et al 2023 .
     */

class ToddlersSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the maximum PAH to Dust fraction value. */
    ENUM_DEF(PAHfraction, High, Low)
        ENUM_VAL(PAHfraction, High, "High PAH fraction")
        ENUM_VAL(PAHfraction, Low, "Low PAH fraction")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution. */
    ENUM_DEF(Resolution, Low, High)
        ENUM_VAL(Resolution, Low, "Low wavelength resolution (Cloudy Default, R=300)")
        ENUM_VAL(Resolution, High, "High wavelength resolution (Lines at R=1e5, Continua at R=300)")
    ENUM_END()

    ITEM_CONCRETE(ToddlersSEDFamily, SEDFamily, "a Toddlers SED family for emission from star-forming regions")

        PROPERTY_ENUM(pahfraction, PAHfraction, "the maximum PAH to Dust fraction value")
        ATTRIBUTE_DEFAULT_VALUE(pahfraction, "High")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")

    ITEM_END()

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit ToddlersSEDFamily(SimulationItem* parent, PAHfraction pahfraction, Resolution resolution);

protected:
    /** This function opens the appropriate resource file (in SKIRT stored table format). */
    void setupSelfBefore() override;

    //====================== Other functions =====================

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each of these objects specifies unit information
        and a human-readable descripton for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family. It retrieves this
        range from the underlying stored table. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
        unit of wavelength) for the %SED with the specified parameters at the specified wavelength,
        or zero if the wavelength is outside of the %SED's intrinsic wavelength range. The number
        and type of parameters must match the information returned by the parameterInfo() function;
        if not the behavior is undefined. */
    double specificLuminosity(double wavelength, const Array& parameters) const override;

    /** This function constructs both the normalized probability density function (pdf) and the
        corresponding normalized cumulative distribution function (cdf) for the %SED with the
        specified parameters over the specified wavelength range. The function returns the
        normalization factor. The number and type of parameters must match the information returned
        by the parameterInfo() function; if not the behavior is undefined. */
    double cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
               const Array& parameters) const override;

    //====================== Data members =====================

private:
    StoredTable<5> _table;
};

////////////////////////////////////////////////////////////////////

#endif

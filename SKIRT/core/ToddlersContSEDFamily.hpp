/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSCONTSEDFAMILY_HPP
#define TODDLERSCONTSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the ToddlersContSEDFamily class represents the family of 
    Toddlers Hii region template SEDs, parameterized on age, metallicity,
    star formation efficiency, and cloud number density and scaled by particle mass.
    
    TODDLERS = Time evolution of Dust Diagnostics and Line Emission from Regions
    containing young Stars.

    The scaling assumes a cloud mass function obeying a power law with slope of -1.8, i.e.,
    dN/dM \propto M^{-1.8} running from cloud mass of 1e5 - 1e6.5.
    The model spectra have 3666 wavelength bins and cover 75 ages from .1-30Myr, 5 metallicity values from
    .001 to .04 (solar being .014), 7 star formation efficiencies from 1 - 15 %, and 8 cloud number
    densities 
    The line and the continuum emission are currently separated, this class is 
    used for the continuum part (stellar + nebular) of the SED.

    The SEDs are tabulated over a wavelength range from 100 Angstrom to 2000 micron.
    The data was generated using the Cloudy master branch (July 2022).
     */

class ToddlersContSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(ToddlersContSEDFamily, SEDFamily, "a Toddlers SED family for continuum emission from Hii regions")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ToddlersContSEDFamily, "Level2")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit ToddlersContSEDFamily(SimulationItem* parent);

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

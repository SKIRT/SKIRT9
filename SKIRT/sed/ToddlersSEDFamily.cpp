/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersSEDFamily.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

ToddlersSEDFamily::ToddlersSEDFamily(SimulationItem* parent, SedMode sedMode, StellarTemplate stellarTemplate,
                                     bool includeDust, Resolution resolution, SFRPeriod sfrPeriod)
{
    parent->addChild(this);
    _sedMode = sedMode;
    _stellarTemplate = stellarTemplate;
    _includeDust = includeDust;
    _resolution = resolution;
    _sfrPeriod = sfrPeriod;
    setup();
}

////////////////////////////////////////////////////////////////////

void ToddlersSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    if (_sedMode == SedMode::SFRNormalizedVariableDust && _stellarTemplate != StellarTemplate::BPASSChab100Bin)
        throw FATALERROR(
            "The SFRNormalizedVariableDust mode is available only for the BPASSChab100Bin stellar template");

    // --- construct the name of the resource file corresponding to the property settings ---

    // base name
    string name = "ToddlersSEDFamily";

    // SED mode
    switch (_sedMode)
    {
        case SedMode::Cloud: name += "_Cloud"; break;
        case SedMode::SFRNormalized: name += "_SFRNormalized"; break;
        case SedMode::SFRNormalizedVariableDust: name += "_SFRNormalizedVariableDust"; break;
    }

    // stellar population
    switch (_stellarTemplate)
    {
        case StellarTemplate::SB99Kroupa100Sin: name += "_SB99_kroupa100_sin"; break;
        case StellarTemplate::BPASSChab100Bin: name += "_BPASS_chab100_bin"; break;
        case StellarTemplate::BPASSChab300Bin: name += "_BPASS_chab300_bin"; break;
    }

    // dust option
    name += _includeDust ? "_Dust" : "_noDust";

    // resolution
    name += _resolution == Resolution::Low ? "_lr" : "_hr";

    // integration period (only for the SFR-normalized modes)
    if (_sedMode == SedMode::SFRNormalized || _sedMode == SedMode::SFRNormalizedVariableDust)
    {
        name += _sfrPeriod == SFRPeriod::Period30Myr ? "_30Myr" : "_10Myr";
    }

    // --- open the resource file ---

    switch (_sedMode)
    {
        case SedMode::Cloud:
            _cloudTable.open(this, name, "lambda(m),time(Myr),Z(1),SFE(1),n_cl(1/cm3),M_cl(Msun)", "Llambda(W/m)",
                             false);
            break;
        case SedMode::SFRNormalized:
            _sfrNormalizedTable.open(this, name, "lambda(m),Z(1),SFE(1),n_cl(1/cm3)", "Llambda(W/m)", false);
            break;
        case SedMode::SFRNormalizedVariableDust:
            _sfrNormalizedVarDustTable.open(this, name, "lambda(m),Z(1),SFE(1),n_cl(1/cm3),f_dust(1)", "Llambda(W/m)",
                                            false);
            break;
    }
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> ToddlersSEDFamily::parameterInfo() const
{
    if (_sedMode == SedMode::Cloud)
    {
        return {
            SnapshotParameter::age(),
            SnapshotParameter::metallicity(),
            SnapshotParameter::custom("star formation efficiency"),
            SnapshotParameter::custom("cloud number density", "numbervolumedensity", "1/cm3"),
            SnapshotParameter::custom("cloud mass", "mass", "Msun"),
            SnapshotParameter::custom("scaling factor"),
        };
    }
    else if (_sedMode == SedMode::SFRNormalizedVariableDust)
    {
        return {
            SnapshotParameter::metallicity(),
            SnapshotParameter::custom("star formation efficiency"),
            SnapshotParameter::custom("cloud number density", "numbervolumedensity", "1/cm3"),
            SnapshotParameter::custom("grain abundance scaling factor"),
            SnapshotParameter::custom("star formation rate", "massrate", "Msun/yr"),
        };
    }
    else  // SFRNormalized
    {
        return {
            SnapshotParameter::metallicity(),
            SnapshotParameter::custom("star formation efficiency"),
            SnapshotParameter::custom("cloud number density", "numbervolumedensity", "1/cm3"),
            SnapshotParameter::custom("star formation rate", "massrate", "Msun/yr"),
        };
    }
}

////////////////////////////////////////////////////////////////////

Range ToddlersSEDFamily::intrinsicWavelengthRange() const
{
    if (_sedMode == SedMode::Cloud)
    {
        return _cloudTable.axisRange<0>();
    }
    else if (_sedMode == SedMode::SFRNormalizedVariableDust)
    {
        return _sfrNormalizedVarDustTable.axisRange<0>();
    }
    else  // SFRNormalized
    {
        return _sfrNormalizedTable.axisRange<0>();
    }
}

////////////////////////////////////////////////////////////////////

double ToddlersSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    if (_sedMode == SedMode::Cloud)
    {
        // Extract parameters for Cloud mode
        double age = parameters[0] / (1e6 * Constants::year());  // Convert from s to Myr
        double Z = parameters[1];                                // Metallicity
        double SFE = parameters[2];                              // Star formation efficiency
        double n_cl = parameters[3] / 1e6;                       // Convert from 1/m³ to 1/cm³
        double M_cl = parameters[4] / Constants::Msun();         // Convert from kg to Msun
        double scaling = parameters[5];                          // Optional scaling factor

        // Use 6D table
        return scaling * _cloudTable(wavelength, age, Z, SFE, n_cl, M_cl);
    }
    else if (_sedMode == SedMode::SFRNormalizedVariableDust)
    {
        // extract parameters for SFRNormalized mode with a grain abundance scaling factor
        double Z = parameters[0];
        double SFE = parameters[1];
        double n_cl = parameters[2] / 1e6;  // Convert from 1/m³ to 1/cm³
        double f_dust = parameters[3];
        double sfr = parameters[4] / Constants::Msun() * Constants::year();  // Convert from kg/s to Msun/yr

        // use 5D table
        return sfr * _sfrNormalizedVarDustTable(wavelength, Z, SFE, n_cl, f_dust);
    }
    else  // SFRNormalized
    {
        // extract parameters for SFRNormalized mode
        double Z = parameters[0];
        double SFE = parameters[1];
        double n_cl = parameters[2] / 1e6;                                   // Convert from 1/m³ to 1/cm³
        double sfr = parameters[3] / Constants::Msun() * Constants::year();  // Convert from kg/s to Msun/yr

        // use 4D table
        return sfr * _sfrNormalizedTable(wavelength, Z, SFE, n_cl);
    }
}

////////////////////////////////////////////////////////////////////

double ToddlersSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                              const Array& parameters) const
{
    if (_sedMode == SedMode::Cloud)
    {
        // Extract parameters for Cloud mode
        double age = parameters[0] / (1e6 * Constants::year());  // Convert from s to Myr
        double Z = parameters[1];                                // Metallicity
        double SFE = parameters[2];                              // Star formation efficiency
        double n_cl = parameters[3] / 1e6;                       // Convert from 1/m³ to 1/cm³
        double M_cl = parameters[4] / Constants::Msun();         // Convert from kg to Msun
        double scaling = parameters[5];                          // Optional scaling factor

        // Use 6D table
        return scaling * _cloudTable.cdf(lambdav, pv, Pv, wavelengthRange, age, Z, SFE, n_cl, M_cl);
    }
    else if (_sedMode == SedMode::SFRNormalizedVariableDust)
    {
        // Extract parameters for SFRNormalized mode with a grain abundance scaling factor
        double Z = parameters[0];
        double SFE = parameters[1];
        double n_cl = parameters[2] / 1e6;  // Convert from 1/m³ to 1/cm³
        double f_dust = parameters[3];
        double sfr = parameters[4] / Constants::Msun() * Constants::year();  // Convert from kg/s to Msun/yr

        // Use 5D table
        return sfr * _sfrNormalizedVarDustTable.cdf(lambdav, pv, Pv, wavelengthRange, Z, SFE, n_cl, f_dust);
    }
    else  // SFRNormalized
    {
        // Extract parameters for SFRNormalized mode
        double Z = parameters[0];
        double SFE = parameters[1];
        double n_cl = parameters[2] / 1e6;                                   // Convert from 1/m³ to 1/cm³
        double sfr = parameters[3] / Constants::Msun() * Constants::year();  // Convert from kg/s to Msun/yr

        // Use 4D table
        return sfr * _sfrNormalizedTable.cdf(lambdav, pv, Pv, wavelengthRange, Z, SFE, n_cl);
    }
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersContSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

ToddlersContSEDFamily::ToddlersContSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void ToddlersContSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "TODDLERS_cont_CMFslope=-1.8", "lambda(m),t(Myr),Z(1),SFE(1),n_cl(1/cm3)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> ToddlersContSEDFamily::parameterInfo() const
{
    return {
        SnapshotParameter::age(),
        SnapshotParameter::metallicity(),
        SnapshotParameter::custom("Star formation efficiency"),
        SnapshotParameter::custom("Cloud number density", "numbervolumedensity", "1/cm3"),
        SnapshotParameter::custom("Stellar mass", "mass", "Msun"),
    };
}

////////////////////////////////////////////////////////////////////

Range ToddlersContSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double ToddlersContSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const

{
    double age     = parameters[0] / (1e6*Constants::year()); 
    double Z       = parameters[1];
    double SFE     = parameters[2];
    double n_cl    = parameters[3] / 1e6 ;
    double M       = parameters[4] / Constants::Msun();
    
    return M * _table(wavelength, age, Z, SFE, n_cl);
}

////////////////////////////////////////////////////////////////////

double ToddlersContSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                              const Array& parameters) const

{
    double age     = parameters[0] / (1e6*Constants::year()); 
    double Z       = parameters[1];
    double SFE     = parameters[2];
    double n_cl    = parameters[3] / 1e6;
    double M       = parameters[4] / Constants::Msun();


    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, age, Z, SFE, n_cl);
}
////////////////////////////////////////////////////////////////////

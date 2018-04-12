/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Units.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

string Units::unit(string qty) const
{
    return _unitDef.unit(qty, type());
}

////////////////////////////////////////////////////////////////////

double Units::in(std::string qty, std::string unit, double value) const
{
    return _unitDef.in(qty, unit, value);
}

////////////////////////////////////////////////////////////////////

double Units::out(string qty, double value) const
{
    return _unitDef.out(qty, type(), value);
}

////////////////////////////////////////////////////////////////////

string Units::ulength() const
{
    return unit("length");
}

////////////////////////////////////////////////////////////////////

double Units::olength(double x) const
{
    return out("length", x);
}

////////////////////////////////////////////////////////////////////

string Units::udistance() const
{
    return unit("distance");
}

////////////////////////////////////////////////////////////////////

double Units::odistance(double d) const
{
    return out("distance", d);
}

////////////////////////////////////////////////////////////////////

string Units::uwavelength() const
{
    return unit("wavelength");
}

////////////////////////////////////////////////////////////////////

double Units::owavelength(double lambda) const
{
    return out("wavelength", lambda);
}

////////////////////////////////////////////////////////////////////

string Units::ugrainsize() const
{
    return unit("grainsize");
}

////////////////////////////////////////////////////////////////////

double Units::ograinsize(double a) const
{
    return out("grainsize", a);
}

////////////////////////////////////////////////////////////////////

string Units::usection() const
{
    return unit("section");
}

////////////////////////////////////////////////////////////////////

double Units::osection(double C) const
{
    return out("section", C);
}

////////////////////////////////////////////////////////////////////

string Units::uvolume() const
{
    return unit("volume");
}

////////////////////////////////////////////////////////////////////

double Units::ovolume(double V) const
{
    return out("volume", V);
}

////////////////////////////////////////////////////////////////////

string Units::uvelocity() const
{
    return unit("velocity");
}

////////////////////////////////////////////////////////////////////

double Units::ovelocity(double v) const
{
    return out("velocity", v);
}

////////////////////////////////////////////////////////////////////

string Units::umass() const
{
    return unit("mass");
}

////////////////////////////////////////////////////////////////////

double Units::omass(double M) const
{
    return out("mass", M);
}

////////////////////////////////////////////////////////////////////

string Units::ubulkmass() const
{
    return unit("bulkmass");
}

////////////////////////////////////////////////////////////////////

double Units::obulkmass(double mu) const
{
    return out("bulkmass", mu);
}

////////////////////////////////////////////////////////////////////

string Units::ubulkmassdensity() const
{
    return unit("bulkmassdensity");
}

////////////////////////////////////////////////////////////////////

double Units::obulkmassdensity(double rho) const
{
    return out("bulkmassdensity", rho);
}

////////////////////////////////////////////////////////////////////

string Units::umasssurfacedensity() const
{
    return unit("masssurfacedensity");
}

////////////////////////////////////////////////////////////////////

double Units::omasssurfacedensity(double Sigma) const
{
    return out("masssurfacedensity", Sigma);
}

////////////////////////////////////////////////////////////////////

string Units::umassvolumedensity() const
{
    return unit("massvolumedensity");
}

////////////////////////////////////////////////////////////////////

double Units::omassvolumedensity(double rho) const
{
    return out("massvolumedensity", rho);
}

////////////////////////////////////////////////////////////////////

string Units::uopacity() const
{
    return unit("opacity");
}

////////////////////////////////////////////////////////////////////

double Units::oopacity(double kappa) const
{
    return out("opacity", kappa);
}

////////////////////////////////////////////////////////////////////

string Units::uenergy() const
{
    return unit("energy");
}

////////////////////////////////////////////////////////////////////

double Units::oenergy(double E) const
{
    return out("energy", E);
}

////////////////////////////////////////////////////////////////////

string Units::ubolluminosity() const
{
    return unit("bolluminosity");
}

////////////////////////////////////////////////////////////////////

double Units::obolluminosity(double L) const
{
    return out("bolluminosity", L);
}

////////////////////////////////////////////////////////////////////

string Units::umonluminosity() const
{
    return unit("monluminosity");
}

////////////////////////////////////////////////////////////////////

double Units::omonluminosity(double Llambda) const
{
    return out("monluminosity", Llambda);
}

////////////////////////////////////////////////////////////////////

string Units::sfluxdensity() const
{
    switch (_fluxOutputStyle)
    {
    case FluxOutputStyle::Wavelength: return "F_lambda";
    case FluxOutputStyle::Frequency:  return "F_nu";
    case FluxOutputStyle::Neutral:    return "lambda*F_lambda";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

string Units::ufluxdensity() const
{
    switch (_fluxOutputStyle)
    {
    case FluxOutputStyle::Wavelength: return unit("wavelengthfluxdensity");
    case FluxOutputStyle::Frequency:  return unit("frequencyfluxdensity");
    case FluxOutputStyle::Neutral:    return unit("neutralfluxdensity");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::ofluxdensity(double lambda, double Flambda) const
{
    switch (_fluxOutputStyle)
    {
    case FluxOutputStyle::Wavelength: return out("wavelengthfluxdensity", Flambda);
    case FluxOutputStyle::Frequency:  return out("frequencyfluxdensity", lambda*lambda*Flambda/Constants::c());
    case FluxOutputStyle::Neutral:    return out("neutralfluxdensity", lambda*Flambda);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

string Units::ssurfacebrightness() const
{
    switch (_fluxOutputStyle)
    {
    case FluxOutputStyle::Wavelength: return "f_lambda";
    case FluxOutputStyle::Frequency:  return "f_nu";
    case FluxOutputStyle::Neutral:    return "lambda*f_lambda";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

string Units::usurfacebrightness() const
{
    switch (_fluxOutputStyle)
    {
    case FluxOutputStyle::Wavelength: return unit("wavelengthsurfacebrightness");
    case FluxOutputStyle::Frequency:  return unit("frequencysurfacebrightness");
    case FluxOutputStyle::Neutral:    return unit("neutralsurfacebrightness");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::osurfacebrightness(double lambda, double flambda) const
{
    switch (_fluxOutputStyle)
    {
    case FluxOutputStyle::Wavelength: return out("wavelengthsurfacebrightness", flambda);
    case FluxOutputStyle::Frequency:  return out("frequencysurfacebrightness", lambda*lambda*flambda/Constants::c());
    case FluxOutputStyle::Neutral:    return out("neutralsurfacebrightness", lambda*flambda);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

string Units::utemperature() const
{
    return unit("temperature");
}

////////////////////////////////////////////////////////////////////

double Units::otemperature(double T) const
{
    return out("temperature", T);
}

////////////////////////////////////////////////////////////////////

string Units::uangle() const
{
    return unit("angle");
}

////////////////////////////////////////////////////////////////////

double Units::oangle(double theta) const
{
    return out("angle", theta);
}

////////////////////////////////////////////////////////////////////

string Units::uposangle() const
{
    return unit("posangle");
}

////////////////////////////////////////////////////////////////////

double Units::oposangle(double theta) const
{
    return out("posangle", theta);
}

////////////////////////////////////////////////////////////////////

string Units::usolidangle() const
{
    return unit("solidangle");
}

////////////////////////////////////////////////////////////////////

double Units::osolidangle(double Omega) const
{
    return out("solidangle", Omega);
}

////////////////////////////////////////////////////////////////////

string Units::upressure() const
{
    return unit("pressure");
}

////////////////////////////////////////////////////////////////////

double Units::opressure(double p) const
{
    return out("pressure", p);
}

////////////////////////////////////////////////////////////////////

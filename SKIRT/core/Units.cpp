/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Units.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    const double c = Constants::c();
    const double hc2 = Constants::h() * Constants::c() * Constants::h() * Constants::c();
}

////////////////////////////////////////////////////////////////////

bool Units::has(string qty, string unit) const
{
    return _unitDef.has(qty, unit);
}

////////////////////////////////////////////////////////////////////

std::tuple<double, double, double> Units::def(string qty, string unit) const
{
    return _unitDef.def(qty, unit);
}

////////////////////////////////////////////////////////////////////

double Units::in(string qty, std::string unit, double value) const
{
    return _unitDef.in(qty, unit, value);
}

////////////////////////////////////////////////////////////////////

double Units::fromFluxStyle(double lambda, double L, FluxOutputStyle style)
{
    switch (style)
    {
        case FluxOutputStyle::Neutral: return L / lambda;
        case FluxOutputStyle::Wavelength: return L;
        case FluxOutputStyle::Frequency: return L * c / lambda / lambda;
        case FluxOutputStyle::Energy: return L * hc2 / lambda / lambda / lambda;
    }
    return L;
}

////////////////////////////////////////////////////////////////////

Array Units::fromFluxStyle(const Array& lambdav, const Array& Lv, FluxOutputStyle style)
{
    switch (style)
    {
        case FluxOutputStyle::Neutral: return Lv / lambdav;
        case FluxOutputStyle::Wavelength: return Lv;
        case FluxOutputStyle::Frequency: return Lv * c / lambdav / lambdav;
        case FluxOutputStyle::Energy: return Lv * hc2 / lambdav / lambdav / lambdav;
    }
    return Lv;
}

////////////////////////////////////////////////////////////////////

string Units::unit(string qty) const
{
    return _unitDef.unit(qty, type());
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

string Units::swavelength() const
{
    switch (_wavelengthOutputStyle)
    {
        case WavelengthOutputStyle::Wavelength: return "lambda";
        case WavelengthOutputStyle::Frequency: return "nu";
        case WavelengthOutputStyle::Energy: return "E";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

bool Units::rwavelength() const
{
    switch (_wavelengthOutputStyle)
    {
        case WavelengthOutputStyle::Wavelength: return false;
        case WavelengthOutputStyle::Frequency: return true;
        case WavelengthOutputStyle::Energy: return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

string Units::uwavelength() const
{
    switch (_wavelengthOutputStyle)
    {
        case WavelengthOutputStyle::Wavelength: return unit("wavelengthwavelength");
        case WavelengthOutputStyle::Frequency: return unit("frequencywavelength");
        case WavelengthOutputStyle::Energy: return unit("energywavelength");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::owavelength(double lambda) const
{
    switch (_wavelengthOutputStyle)
    {
        case WavelengthOutputStyle::Wavelength: return out("wavelengthwavelength", lambda);
        case WavelengthOutputStyle::Frequency: return out("frequencywavelength", lambda);
        case WavelengthOutputStyle::Energy: return out("energywavelength", lambda);
    }
    return 0.;
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

std::string Units::upergrainsize() const
{
    return unit("pergrainsize");
}

////////////////////////////////////////////////////////////////////

double Units::opergrainsize(double a) const
{
    return out("pergrainsize", a);
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

string Units::unumbersurfacedensity() const
{
    return unit("numbersurfacedensity");
}

////////////////////////////////////////////////////////////////////

double Units::onumbersurfacedensity(double N) const
{
    return out("numbersurfacedensity", N);
}

////////////////////////////////////////////////////////////////////

string Units::unumbervolumedensity() const
{
    return unit("numbervolumedensity");
}

////////////////////////////////////////////////////////////////////

double Units::onumbervolumedensity(double n) const
{
    return out("numbervolumedensity", n);
}

////////////////////////////////////////////////////////////////////

string Units::umasscoefficient() const
{
    return unit("masscoefficient");
}

////////////////////////////////////////////////////////////////////

double Units::omasscoefficient(double kappa) const
{
    return out("masscoefficient", kappa);
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

std::string Units::umagneticfield() const
{
    return unit("magneticfield");
}

////////////////////////////////////////////////////////////////////

double Units::omagneticfield(double B) const
{
    return out("magneticfield", B);
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

string Units::ubolluminosityvolumedensity() const
{
    return unit("bolluminosityvolumedensity");
}

////////////////////////////////////////////////////////////////////

double Units::obolluminosityvolumedensity(double L) const
{
    return out("bolluminosityvolumedensity", L);
}

////////////////////////////////////////////////////////////////////

string Units::ubolluminositysurfacedensity() const
{
    return unit("bolluminositysurfacedensity");
}

////////////////////////////////////////////////////////////////////

double Units::obolluminositysurfacedensity(double L) const
{
    return out("bolluminositysurfacedensity", L);
}

////////////////////////////////////////////////////////////////////

string Units::smonluminosity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return "lambda*L_lambda";
        case FluxOutputStyle::Wavelength: return "L_lambda";
        case FluxOutputStyle::Frequency: return "L_nu";
        case FluxOutputStyle::Energy: return "L_E";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

string Units::umonluminosity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return unit("neutralmonluminosity");
        case FluxOutputStyle::Wavelength: return unit("wavelengthmonluminosity");
        case FluxOutputStyle::Frequency: return unit("frequencymonluminosity");
        case FluxOutputStyle::Energy: return unit("energymonluminosity");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::omonluminosity(double lambda, double Llambda) const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return out("neutralmonluminosity", lambda * Llambda);
        case FluxOutputStyle::Wavelength: return out("wavelengthmonluminosity", Llambda);
        case FluxOutputStyle::Frequency: return out("frequencymonluminosity", lambda * lambda * Llambda / c);
        case FluxOutputStyle::Energy: return out("energymonluminosity", lambda * lambda * lambda * Llambda / hc2);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

string Units::smonluminosityvolumedensity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return "lambda*L_lambda/V";
        case FluxOutputStyle::Wavelength: return "L_lambda/V";
        case FluxOutputStyle::Frequency: return "L_nu/V";
        case FluxOutputStyle::Energy: return "L_E/V";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

string Units::umonluminosityvolumedensity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return unit("neutralmonluminosityvolumedensity");
        case FluxOutputStyle::Wavelength: return unit("wavelengthmonluminosityvolumedensity");
        case FluxOutputStyle::Frequency: return unit("frequencymonluminosityvolumedensity");
        case FluxOutputStyle::Energy: return unit("energymonluminosityvolumedensity");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::omonluminosityvolumedensity(double lambda, double Llambda) const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return out("neutralmonluminosityvolumedensity", lambda * Llambda);
        case FluxOutputStyle::Wavelength: return out("wavelengthmonluminosityvolumedensity", Llambda);
        case FluxOutputStyle::Frequency:
            return out("frequencymonluminosityvolumedensity", lambda * lambda * Llambda / c);
        case FluxOutputStyle::Energy:
            return out("energymonluminosityvolumedensity", lambda * lambda * lambda * Llambda / hc2);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

string Units::sfluxdensity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return "lambda*F_lambda";
        case FluxOutputStyle::Wavelength: return "F_lambda";
        case FluxOutputStyle::Frequency: return "F_nu";
        case FluxOutputStyle::Energy: return "F_E";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

string Units::ufluxdensity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return unit("neutralfluxdensity");
        case FluxOutputStyle::Wavelength: return unit("wavelengthfluxdensity");
        case FluxOutputStyle::Frequency: return unit("frequencyfluxdensity");
        case FluxOutputStyle::Energy: return unit("energyfluxdensity");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::ofluxdensity(double lambda, double Flambda) const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return out("neutralfluxdensity", lambda * Flambda);
        case FluxOutputStyle::Wavelength: return out("wavelengthfluxdensity", Flambda);
        case FluxOutputStyle::Frequency: return out("frequencyfluxdensity", lambda * lambda * Flambda / c);
        case FluxOutputStyle::Energy: return out("energyfluxdensity", lambda * lambda * lambda * Flambda / hc2);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

string Units::ssurfacebrightness() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return "lambda*f_lambda";
        case FluxOutputStyle::Wavelength: return "f_lambda";
        case FluxOutputStyle::Frequency: return "f_nu";
        case FluxOutputStyle::Energy: return "f_E";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

string Units::usurfacebrightness() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return unit("neutralsurfacebrightness");
        case FluxOutputStyle::Wavelength: return unit("wavelengthsurfacebrightness");
        case FluxOutputStyle::Frequency: return unit("frequencysurfacebrightness");
        case FluxOutputStyle::Energy: return unit("energysurfacebrightness");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::osurfacebrightness(double lambda, double flambda) const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return out("neutralsurfacebrightness", lambda * flambda);
        case FluxOutputStyle::Wavelength: return out("wavelengthsurfacebrightness", flambda);
        case FluxOutputStyle::Frequency: return out("frequencysurfacebrightness", lambda * lambda * flambda / c);
        case FluxOutputStyle::Energy: return out("energysurfacebrightness", lambda * lambda * lambda * flambda / hc2);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

std::string Units::smeanintensity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return "lambda*J_lambda";
        case FluxOutputStyle::Wavelength: return "J_lambda";
        case FluxOutputStyle::Frequency: return "J_nu";
        case FluxOutputStyle::Energy: return "J_E";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

std::string Units::umeanintensity() const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return unit("neutralmeanintensity");
        case FluxOutputStyle::Wavelength: return unit("wavelengthmeanintensity");
        case FluxOutputStyle::Frequency: return unit("frequencymeanintensity");
        case FluxOutputStyle::Energy: return unit("energymeanintensity");
    }
    return string();
}

////////////////////////////////////////////////////////////////////

double Units::omeanintensity(double lambda, double Jlambda) const
{
    switch (_fluxOutputStyle)
    {
        case FluxOutputStyle::Neutral: return out("neutralmeanintensity", lambda * Jlambda);
        case FluxOutputStyle::Wavelength: return out("wavelengthmeanintensity", Jlambda);
        case FluxOutputStyle::Frequency: return out("frequencymeanintensity", lambda * lambda * Jlambda / c);
        case FluxOutputStyle::Energy: return out("energymeanintensity", lambda * lambda * lambda * Jlambda / hc2);
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

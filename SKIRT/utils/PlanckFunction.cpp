/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PlanckFunction.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

PlanckFunction::PlanckFunction(double T)
    : _T(T)
{
}

////////////////////////////////////////////////////////////////////

double PlanckFunction::operator()(double lambda) const
{
    const double& h = Constants::h();
    const double& c = Constants::c();
    const double& k = Constants::k();
    double x = h*c/(lambda*k*_T);
    return 2.0*h*c*c / pow(lambda,5) / (exp(x)-1.0);
}

////////////////////////////////////////////////////////////////////

double PlanckFunction::derivative(double lambda) const
{
    const double& h = Constants::h();
    const double& c = Constants::c();
    const double& k = Constants::k();
    double x = h*c/(lambda*k*_T);
    double ans = (2.0*h*c*c*x/_T) / pow(lambda,5) / (2.0*(cosh(x)-1.0));
    return ans;
}

////////////////////////////////////////////////////////////////////

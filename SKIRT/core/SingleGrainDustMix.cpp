/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SingleGrainDustMix.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StokesVector.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // arbitrary constants determining the resolution for disretizing the scattering angles
    constexpr int Ntheta = 181;                 // theta from 0 to pi, index t
    constexpr double dtheta = M_PI/(Ntheta-1);
    constexpr int Nphi = 361;                   // phi from 0 to 2 pi, index f
    constexpr double dphi = 2*M_PI/(Nphi-1);
}

////////////////////////////////////////////////////////////////////

void SingleGrainDustMix::setupSelfBefore()
{
    MaterialMix::setupSelfBefore();

    // get the basic optical properties
    string name = resourceNameForOpticalProps();
    _sigmaabs.open(this, name, "lambda(m)", "sigmaabs(m2/H)");
    _sigmasca.open(this, name, "lambda(m)", "sigmasca(m2/H)");
    _asymmpar.open(this, name, "lambda(m)", "g(1)");
    _mu = StoredTable<1>(this, name, "lambda(m)", "mu(kg/H)")[1.];  // get mu value for arbitrary wavelength

    // get the Mueller matrix elements, as required by the scattering mode
    name = resourceNameForMuellerMatrix();
    if (!name.empty())
    {
        if (scatteringMode() == ScatteringMode::MaterialPhaseFunction ||
            scatteringMode() == ScatteringMode::SphericalPolarization)
        {
            // open the resource for the first Mueller matrix coefficient, which describes the phase function
            _S11.open(this, name, "lambda(m),theta(rad)", "S11(1)");

            // create tables listing theta and sin(theta)*dtheta for a number of t indices
            _thetav.resize(Ntheta);
            _thetasv.resize(Ntheta);
            for (int t=0; t<Ntheta; t++)
            {
                double theta = t * dtheta;
                _thetav[t] = theta;
                _thetasv[t] = sin(theta) * dtheta;
            }

            if (scatteringMode() == ScatteringMode::SphericalPolarization)
            {
                // open the resources for the other Mueller matrix coefficients as well
                _S12.open(this, name, "lambda(m),theta(rad)", "S12(1)");
                _S33.open(this, name, "lambda(m),theta(rad)", "S33(1)");
                _S34.open(this, name, "lambda(m),theta(rad)", "S34(1)");

                // create tables listing phi, phi/(2 pi), sin(2 phi) and 1-cos(2 phi) for a number of f indices
                _phiv.resize(Nphi);
                _phi1v.resize(Nphi);
                _phisv.resize(Nphi);
                _phicv.resize(Nphi);
                for (int f=0; f<Nphi; f++)
                {
                    double phi = f * dphi;
                    _phiv[f] = phi;
                    _phi1v[f] = phi/(2*M_PI);
                    _phisv[f] = sin(2*phi);
                    _phicv[f] = 1-cos(2*phi);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

string SingleGrainDustMix::resourceNameForMuellerMatrix() const
{
    return string();
}

////////////////////////////////////////////////////////////////////

MaterialMix::MaterialType SingleGrainDustMix::materialType() const
{
    return MaterialType::Dust;
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::mass() const
{
    return _mu;
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::sectionAbsSelf(double lambda) const
{
    return _sigmaabs[lambda];
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::sectionScaSelf(double lambda) const
{
    return _sigmasca[lambda];
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::asymmpar(double lambda) const
{
    return _asymmpar[lambda];
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::phaseFunctionValueForCosine(double lambda, double costheta) const
{
    // calculate the phase function normalization factor for this wavelength
    double sum = 0.;        // (sin(0) and sin(pi) are zero anyway)
    for (int t=1; t<Ntheta-1; t++) sum += _S11(lambda,_thetav[t])*_thetasv[t];
    double normalization = 2.0/sum;

    // calculate the phase function value
    return normalization * _S11(lambda,acos(costheta));
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::generateCosineFromPhaseFunction(double lambda) const
{
    // construct and sample from the normalized cumulative distribution of theta for this wavelength
    Array thetaXv;
    NR::cdf(thetaXv, Ntheta-1, [this,lambda](int t){ return _S11(lambda,_thetav[t+1])*_thetasv[t+1]; });
    return cos(random()->cdfLinLin(_thetav, thetaXv));
}

////////////////////////////////////////////////////////////////////

double SingleGrainDustMix::phaseFunctionValue(double lambda, double theta, double phi, const StokesVector* sv) const
{
    // calculate the phase function normalization factor for this wavelength
    double sum = 0.;        // (sin(0) and sin(pi) are zero anyway)
    for (int t=1; t<Ntheta-1; t++) sum += _S11(lambda,_thetav[t])*_thetasv[t];
    double normalization = 2.0/sum;

    // calculate the phase function value
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    return normalization * (_S11(lambda,theta) + polDegree*_S12(lambda,theta)*cos(2.*(phi-polAngle)));
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> SingleGrainDustMix::generateAnglesFromPhaseFunction(double lambda, const StokesVector* sv) const
{
    // construct and sample from the normalized cumulative distribution of theta for this wavelength
    Array thetaXv;
    NR::cdf(thetaXv, Ntheta-1, [this,lambda](int t){ return _S11(lambda,_thetav[t+1])*_thetasv[t+1]; });
    double theta = random()->cdfLinLin(_thetav, thetaXv);

    // construct and sample from the normalized cumulative distribution of phi for this wavelength and theta angle
    double polDegree = sv->linearPolarizationDegree();
    double polAngle = sv->polarizationAngle();
    double PF = polDegree * _S12(lambda,theta)/_S11(lambda,theta) / (4*M_PI);
    double cos2polAngle = cos(2*polAngle) * PF;
    double sin2polAngle = sin(2*polAngle) * PF;
    double phi = random()->cdfLinLin(_phiv, _phi1v + cos2polAngle*_phisv + sin2polAngle*_phicv);

    // return the result
    return std::make_pair(theta, phi);
}

////////////////////////////////////////////////////////////////////

void SingleGrainDustMix::applyMueller(double lambda, double theta, StokesVector* sv) const
{
    sv->applyMueller(_S11(lambda,theta), _S12(lambda,theta), _S33(lambda,theta), _S34(lambda,theta));
}

////////////////////////////////////////////////////////////////////

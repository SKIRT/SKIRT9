/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiGaussianExpansionGeometry.hpp"
#include "FatalError.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"

//////////////////////////////////////////////////////////////////////

void MultiGaussianExpansionGeometry::setupSelfBefore()
{
    AxGeometry::setupSelfBefore();

    // read the file with the raw MGE parameters
    TextInFile infile(this, _filename, "multi-gaussian expansion parameters");
    infile.addColumn("count or weight");
    infile.addColumn("scalelength in pixels");
    infile.addColumn("apparent flattening");
    infile.readAllColumns(_Mv, _sigmav, _qv);
    infile.close();
    _Ncomp = _Mv.size();
    if (!_Ncomp) throw FATALERROR("Multi-gaussian expansion must have at least one component");

    // convert from pixelscale to physical scale
    _sigmav *= _pixelscale;

    // convert the apparent flattening to real flattening (see e.g. Bacon 1985, A&A, 143, 84)
    double cosi = cos(_inclination);
    double sini = sin(_inclination);
    for (int i = 0; i < _Ncomp; i++)
    {
        if (_qv[i] < cosi)
            throw FATALERROR("MGE component with index " + std::to_string(i)
                             + " can't be deprojected: apparent flattening is smaller than cosine of inclination ("
                             + StringUtils::toString(_qv[i], 'f') + " < " + StringUtils::toString(cosi, 'f') + ")");
        _qv[i] = sqrt((_qv[i] - cosi) * (_qv[i] + cosi)) / sini;
    }

    // convert the counts to normalized weight and set up a vector with cumulative weights
    double Mtot = NR::cdf(_Mcumv, _Mv);
    _Mv /= Mtot;
}

////////////////////////////////////////////////////////////////////

double MultiGaussianExpansionGeometry::density(double R, double z) const
{
    double rho = 0.0;
    for (int i = 0; i < _Ncomp; i++)
    {
        double q = _qv[i];
        double M = _Mv[i];
        double sigma = _sigmav[i];
        double rho0 = M / pow(sqrt(2.0 * M_PI) * sigma, 3) / q;
        double m2 = R * R + z * z / (q * q);
        double sigma2 = sigma * sigma;
        rho += rho0 * exp(-0.5 * m2 / sigma2);
    }
    return rho;
}

////////////////////////////////////////////////////////////////////

Position MultiGaussianExpansionGeometry::generatePosition() const
{
    // select a component according to its relative contribution
    int i = NR::locateClip(_Mcumv, random()->uniform());

    // generate a position for that component
    double q = _qv[i];
    double sigma = _sigmav[i];
    double x = sigma * random()->gauss();
    double y = sigma * random()->gauss();
    double z = q * sigma * random()->gauss();
    return Position(x, y, z);
}

//////////////////////////////////////////////////////////////////////

double MultiGaussianExpansionGeometry::SigmaR() const
{
    double sum = 0.0;
    for (int i = 0; i < _Ncomp; i++)
    {
        double sigma = _sigmav[i];
        double sigma2 = sigma * sigma;
        sum += _Mv[i] / (4.0 * M_PI) / sigma2 / _qv[i];
    }
    return sum;
}

//////////////////////////////////////////////////////////////////////

double MultiGaussianExpansionGeometry::SigmaZ() const
{
    double sum = 0.0;
    for (int i = 0; i < _Ncomp; i++)
    {
        double sigma = _sigmav[i];
        double sigma2 = sigma * sigma;
        sum += _Mv[i] / (2.0 * M_PI) / sigma2;
    }
    return sum;
}

////////////////////////////////////////////////////////////////////

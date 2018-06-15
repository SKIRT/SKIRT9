/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GeometricSource.hpp"
#include "Constants.hpp"
#include "PhotonPacket.hpp"
#include "Random.hpp"
#include "RedshiftInterface.hpp"
#include "WavelengthRangeInterface.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // An instance of this class offers the redshift interface for the bulk velocity specified in the constructor
    class BulkVelocity : public RedshiftInterface
    {
    private:
        Vec _bfv;
    public:
        BulkVelocity(Vec bfv) : _bfv(bfv) { }
        double redshiftForDirection(Direction bfk) const override { return -Vec::dot(_bfv, bfk) / Constants::c(); }
    };
}

//////////////////////////////////////////////////////////////////////

void GeometricSource::setupSelfBefore()
{
    Source::setupSelfBefore();

    if (velocityX() || velocityY() || velocityZ())
        _bulkvelocity = new BulkVelocity(Vec(velocityX(), velocityY(), velocityZ()));
}

//////////////////////////////////////////////////////////////////////

GeometricSource::~GeometricSource()
{
    delete _bulkvelocity;
}

//////////////////////////////////////////////////////////////////////

int GeometricSource::dimension() const
{
    int velocityDimension = 1;
    if (velocityZ()) velocityDimension = 2;
    if (velocityX() || velocityY()) velocityDimension = 3;
    return max(_geometry->dimension(), velocityDimension);
}

//////////////////////////////////////////////////////////////////////

double GeometricSource::luminosity() const
{
    return _normalization->luminosity(_sed);
}

//////////////////////////////////////////////////////////////////////

double GeometricSource::specificLuminosity(double wavelength) const
{
    if (!interface<WavelengthRangeInterface>()->wavelengthRange().contains(wavelength)) return 0.;
    return _sed->specificLuminosity(wavelength) * luminosity();
}

//////////////////////////////////////////////////////////////////////

void GeometricSource::launch(PhotonPacket* pp, size_t historyIndex, double L) const
{
    // generate a random position from the geometry
    Position bfr = _geometry->generatePosition();

    // generate a random wavelength from the SED and/or from the bias distribution
    double lambda, w;
    double xi = wavelengthBias();
    if (!xi)
    {
        // no biasing -- simply use the intrinsic spectral distribution
        lambda = _sed->generateWavelength();
        w = 1.;
    }
    else
    {
        // biasing -- use one or the other distribution
        if (random()->uniform() > xi) lambda = sed()->generateWavelength();
        else lambda = wavelengthBiasDistribution()->generateWavelength();

        // calculate the compensating weight factor
        double s = sed()->specificLuminosity(lambda);
        if (!s)
        {
            // if the wavelength can't occur in the intrinsic distribution,
            // the weight factor is zero regardless of the probability in the bias distribution
            // (handling this separately also avoids NaNs in the pathological case s=d=0)
            w = 0.;
        }
        else
        {
            // regular composite bias weight
            double b = wavelengthBiasDistribution()->probability(lambda);
            w = s / ((1-xi)*s + xi*b);
        }
    }

    // launch the photon packet with isotropic direction
    pp->launch(historyIndex, lambda, L*w, bfr, random()->direction(), _bulkvelocity);
}

//////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Random.hpp"
#include "Box.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "Position.hpp"

//////////////////////////////////////////////////////////////////////

void Random::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

}

//////////////////////////////////////////////////////////////////////

void Random::initialize()
{
}

//////////////////////////////////////////////////////////////////////

void Random::randomize()
{
}

//////////////////////////////////////////////////////////////////////

double Random::uniform()
{
    return 0;
}

//////////////////////////////////////////////////////////////////////

double Random::cdf(const Array& xv, const Array& Xv)
{
    double X = uniform();
    int i = NR::locateClip(Xv, X);
    return NR::interpolateLinLin(X, Xv[i], Xv[i+1], xv[i], xv[i+1]);
}

//////////////////////////////////////////////////////////////////////

double Random::gauss()
{
    double rsq, v1, v2;
    do {
        v1 = 2.0*uniform()-1.0;
        v2 = 2.0*uniform()-1.0;
        rsq = v1*v1+v2*v2;
    } while (rsq>=1.0 || rsq==0.0);
    return v2*sqrt(-2.0*log(rsq)/rsq);
}

//////////////////////////////////////////////////////////////////////

double Random::expon()
{
    return -log(1.0-uniform());
}

//////////////////////////////////////////////////////////////////////

double Random::exponCutoff(double xmax)
{
    if (xmax==0.0)
        return 0.0;
    else if (xmax<1e-10)
        return uniform()*xmax;
    double x = -log(1.0-uniform()*(1.0-exp(-xmax)));
    while (x>xmax)
    {
        x = -log(1.0-uniform()*(1.0-exp(-xmax)));
    }
    return x;
}

//////////////////////////////////////////////////////////////////////

Direction Random::direction()
{
    double theta = acos(2.0*uniform()-1.0);
    double phi = 2.0*M_PI*uniform();
    return Direction(theta,phi);
}

//////////////////////////////////////////////////////////////////////

Direction Random::direction(Direction bfk, double costheta)
{
    // generate random phi and get the sine and cosine for both angles
    double phi = 2.0*M_PI * uniform();
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    double sintheta = sqrt(fabs((1.0-costheta)*(1.0+costheta)));

    // get the old direction
    double kx, ky, kz;
    bfk.cartesian(kx,ky,kz);

    // determine the new direction
    double kxnew, kynew, kznew;
    if (kz>0.99999)
    {
        kxnew = cosphi * sintheta;
        kynew = sinphi * sintheta;
        kznew = costheta;
    }
    else if (kz<-0.99999)
    {
        kxnew = cosphi * sintheta;
        kynew = sinphi * sintheta;
        kznew = -costheta;
    }
    else
    {
        double root = sqrt((1.0-kz)*(1.0+kz));
        kxnew = sintheta/root*(-kx*kz*cosphi+ky*sinphi) + kx*costheta;
        kynew = -sintheta/root*(ky*kz*cosphi+kx*sinphi) + ky*costheta;
        kznew = root*sintheta*cosphi + kz*costheta;
    }
    return Direction(kxnew,kynew,kznew);
}

//////////////////////////////////////////////////////////////////////

Position Random::position(const Box& box)
{
    // generate the random numbers in separate statements to guarantee evaluation order
    // (function arguments are evaluated in different order depending on the compiler)
    double x = uniform();
    double y = uniform();
    double z = uniform();
    return Position(box.fracPos(x,y,z));
}

//////////////////////////////////////////////////////////////////////

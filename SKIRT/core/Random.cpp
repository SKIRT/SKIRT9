/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Random.hpp"
#include "Box.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "Position.hpp"

//////////////////////////////////////////////////////////////////////

void Random::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _parfac = find<ParallelFactory>();
    int Nthreads = _parfac->maxThreadCount();
    _mtv.resize(Nthreads);
    _mtiv.resize(Nthreads);

    initialize(Nthreads);
}

//////////////////////////////////////////////////////////////////////

void Random::initialize(int Nthreads)
{
    unsigned long seed = _seed;
    for (int thread=0; thread<Nthreads; thread++)
    {
        find<Log>()->info("Initializing random number generator for thread number "
                          + std::to_string(thread) + " with seed " +std::to_string(seed) + "... ");
        vector<unsigned long>& mt = _mtv[thread];
        int& mti = _mtiv[thread];
        mt.resize(624);
        mt[0] = seed & 0xffffffff;
        for (mti=1; mti<624; mti++)
            mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
        ++seed;
    }
}

//////////////////////////////////////////////////////////////////////

void Random::randomize()
{
    find<Log>()->info("Setting different seeds for each process.");
    PeerToPeerCommunicator* comm = find<PeerToPeerCommunicator>();

    int Nthreads = _parfac->maxThreadCount();
    _mtv.resize(Nthreads);      // Because the number of threads can be different during and after the setup
    _mtiv.resize(Nthreads);     // of the simulation.

    _seed = _seed + Nthreads * comm->rank();

    initialize(Nthreads);
}

//////////////////////////////////////////////////////////////////////

double Random::uniform()
{
    int thread = _parfac->currentThreadIndex();
    vector<unsigned long>& mt = _mtv[thread];
    int& mti = _mtiv[thread];
    double ans = 0.0;
    do
    {
        unsigned long y;
        static unsigned long mag01[2]={0x0,0x9908b0df};
        if (mti >= 624)
        {
            int kk;
            for (kk=0;kk<227;kk++)
            {
                y = (mt[kk]&0x80000000)|(mt[kk+1]&0x7fffffff);
                mt[kk] = mt[kk+397] ^ (y >> 1) ^ mag01[y & 0x1];
            }
            for (;kk<624-1;kk++)
            {
                y = (mt[kk]&0x80000000)|(mt[kk+1]&0x7fffffff);
                mt[kk] = mt[kk-227] ^ (y >> 1) ^ mag01[y & 0x1];
            }
            y = (mt[623]&0x80000000)|(mt[0]&0x7fffffff);
            mt[623] = mt[396] ^ (y >> 1) ^ mag01[y & 0x1];
            mti = 0;
        }
        y = mt[mti++];
        y ^= (y>>11);
        y ^= (y<<7) & 0x9d2c5680;
        y ^= (y<<15) & 0xefc60000;
        y ^= (y>>18);
        ans = static_cast<double>(y) / static_cast<unsigned long>(0xffffffff);
    }
    while (ans<=0.0 || ans>=1.0);
    return ans;
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

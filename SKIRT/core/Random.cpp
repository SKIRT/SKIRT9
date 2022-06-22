/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Random.hpp"
#include "Box.hpp"
#include "NR.hpp"
#include "Position.hpp"
#include "SpecialFunctions.hpp"
#include <random>
#include <stack>

//////////////////////////////////////////////////////////////////////

namespace
{
    // This helper class represents a pseudo-random generator. An instance is always constructed as
    // an arbitrary generator, but it can be turned into a predictable generator through setState().
    class Rand
    {
    private:
        // use 64-bit Mersenne twister
        std::mt19937_64 _generator;
        // explicitly exclude zero from range; one is excluded automatically
        std::uniform_real_distribution<double> _distribution{
            std::nextafter(static_cast<double>(0.), static_cast<double>(1.)), 1.};

    public:
        // construct arbitrary generator, seeded with a truly random sequence
        Rand()
        {
            std::random_device r;
            std::seed_seq seedseq{r(), r(), r(), r(), r(), r(), r(), r()};
            _generator.seed(seedseq);
        }

        // turn into predictable generator, seeded with fixed sequence depending on given seed
        void setState(int seed)
        {
            std::seed_seq seedseq{979364188u + seed, 871244425u + seed, 1693909487u + seed, 1290454318u + seed,
                                  210509498u + seed, 542237529u + seed, 3429911442u + seed, 3321294726u + seed};
            _generator.seed(seedseq);
        }

        // get uniform deviate
        double get() { return _distribution(_generator); }
    };

    // allocate a random generator for each thread, constructed when the thread is created
    thread_local Rand _rng;

    // allocate a random generator stack for each thread, constructed when the thread is created;
    // this stack is used solely by the push() and pop() functions
    thread_local std::stack<Rand> _stack;
}

//////////////////////////////////////////////////////////////////////

void Random::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // initialize the generator for this thread
    _rng.setState(seed());
}

//////////////////////////////////////////////////////////////////////

double Random::uniform()
{
    return _rng.get();
}

//////////////////////////////////////////////////////////////////////

double Random::gauss()
{
    double rsq, v1, v2;
    do
    {
        v1 = 2.0 * uniform() - 1.0;
        v2 = 2.0 * uniform() - 1.0;
        rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    return v2 * sqrt(-2.0 * log(rsq) / rsq);
}

//////////////////////////////////////////////////////////////////////

double Random::expon()
{
    return -log(uniform());
}

//////////////////////////////////////////////////////////////////////

double Random::exponCutoff(double xmax)
{
    if (xmax == 0.0)
        return 0.0;
    else if (xmax < 1e-10)
        return uniform() * xmax;
    double x = -log(1.0 - uniform() * (1.0 - exp(-xmax)));
    while (x > xmax)
    {
        x = -log(1.0 - uniform() * (1.0 - exp(-xmax)));
    }
    return x;
}

//////////////////////////////////////////////////////////////////////

Direction Random::direction()
{
    double theta = acos(2.0 * uniform() - 1.0);
    double phi = 2.0 * M_PI * uniform();
    return Direction(theta, phi);
}

//////////////////////////////////////////////////////////////////////

Direction Random::direction(Direction bfk, double costheta)
{
    // generate random phi and get the sine and cosine for both angles
    double phi = 2.0 * M_PI * uniform();
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    double sintheta = sqrt(fabs((1.0 - costheta) * (1.0 + costheta)));

    // get the old direction
    double kx, ky, kz;
    bfk.cartesian(kx, ky, kz);

    // determine the new direction
    double kxnew, kynew, kznew;
    if (kz > 0.99999)
    {
        kxnew = cosphi * sintheta;
        kynew = sinphi * sintheta;
        kznew = costheta;
    }
    else if (kz < -0.99999)
    {
        kxnew = cosphi * sintheta;
        kynew = sinphi * sintheta;
        kznew = -costheta;
    }
    else
    {
        double root = sqrt((1.0 - kz) * (1.0 + kz));
        kxnew = sintheta / root * (-kx * kz * cosphi + ky * sinphi) + kx * costheta;
        kynew = -sintheta / root * (ky * kz * cosphi + kx * sinphi) + ky * costheta;
        kznew = root * sintheta * cosphi + kz * costheta;
    }
    return Direction(kxnew, kynew, kznew);
}

//////////////////////////////////////////////////////////////////////

Position Random::position(const Box& box)
{
    // generate the random numbers in separate statements to guarantee evaluation order
    // (function arguments are evaluated in different order depending on the compiler)
    double x = uniform();
    double y = uniform();
    double z = uniform();
    return Position(box.fracPos(x, y, z));
}

Vec Random::maxwell()
{
    // generate the random numbers in separate statements to guarantee evaluation order
    // (function arguments are evaluated in different order depending on the compiler)
    double x = gauss();
    double y = gauss();
    double z = gauss();
    return Vec(x, y, z);
}

//////////////////////////////////////////////////////////////////////

double Random::cdfLinLin(const Array& xv, const Array& Pv)
{
    double X = uniform();
    int i = NR::locateClip(Pv, X);
    return NR::interpolateLinLin(X, Pv[i], Pv[i + 1], xv[i], xv[i + 1]);
}

//////////////////////////////////////////////////////////////////////

double Random::cdfLogLog(const Array& xv, const Array& pv, const Array& Pv)
{
    double X = uniform();
    int i = NR::locateClip(Pv, X);
    double alpha = log(pv[i + 1] / pv[i]) / log(xv[i + 1] / xv[i]);
    return xv[i] * SpecialFunctions::gexp(-alpha, (X - Pv[i]) / (pv[i] * xv[i]));
}

//////////////////////////////////////////////////////////////////////

void Random::push(int seed)
{
    _stack.push(_rng);
    _rng.setState(seed);
}

//////////////////////////////////////////////////////////////////////

void Random::pop()
{
    _rng = _stack.top();
    _stack.pop();
}

//////////////////////////////////////////////////////////////////////

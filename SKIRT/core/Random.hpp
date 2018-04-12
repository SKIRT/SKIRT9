/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
class Box;
class Direction;
class ParallelFactory;
class Position;

//////////////////////////////////////////////////////////////////////

/** This class contains a random number generator, and can be used to produce series of random
    numbers for different probability distributions. Typically, only a single instance of the class
    should be constructed for each simulation. This random number class is adapted from a C library
    known as genrand(), written by Takuji Nishimura. More information can be found at
    http://www.math.keio.ac.jp/matumoto/emt.html. */
class Random : public SimulationItem
{
    ITEM_CONCRETE(Random, SimulationItem, "the default random generator")

    PROPERTY_INT(seed, "the seed for the random generator")
        ATTRIBUTE_MIN_VALUE(seed, "1")
        ATTRIBUTE_DEFAULT_VALUE(seed, "4357")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function initializes the data structure that is used for the calculation of the random
        numbers, using the current value of the seed. */
    void setupSelfBefore() override;

private:
    /** This function serves as the body of two different functions: the setupSelfBefore() function and
        the randomize() function. Based on the seed stored in the \c _seed attribute, this function
        generates random sequences for each different thread in the simulation. This is done by
        incrementing the seed with each increment of the thread number. During setup, the \c _seed
        variable is equal on each process, providing them with the same random sequences. When this
        function is called from randomize(), each process has given a different seed, yielding
        different random sequences for every thread in the multiprocessing environment. */
    void initialize(int Nthreads);

    //======================== Other Functions =======================

public:
    /** This function is used to give each thread in the multiprocessing environment a different random
        seed and regenerating their random sequences. This function uses the find algorithm to obtain a
        pointer to the PeerToPeerCommunicator object. Next, the \c _seed variable is shifted by exactly
        the number of threads for each successive process. At least, the initialize() function is
        called with the number of threads as an argument. */
    void randomize();

    /** This function generates a random uniform deviate, i.e. a random double precision number in
        the interval [0,1]. For details how this is exactly done, see the information at
        http://www.math.keio.ac.jp/matumoto/emt.html. */
    double uniform();

    /** This function generates a random number drawn from an arbitrary probability distribution
        \f$p(x)\,{\text{d}}x\f$ with corresponding cumulative distribution function \f$P(x)\f$.
        The routine reads in a discretized version \f$P_i\f$ of the cdf sampled at a set of
        points \f$x_i\f$. A uniform deviate \f$X\f$ is generated, and the equation \f$X=P(x)\f$
        is solved using linear interpolation. */
    double cdf(const Array& xv, const Array& Xv);

    /** This function generates a random number from a Gaussian distribution function with mean 0
        and standard deviation 1, i.e. defined by the probability distribution \f[ p(x)\,{\rm d}x =
        \frac{1}{\sqrt{2\pi}}\, {\rm e}^{-\frac12\,x^2}\,{\rm d}x.\f] The algorithm used and the
        implementation are taken from Press et al. (2002). */
    double gauss();

    /** This function generates a random number from an exponential distribution function, defined
        by the probability distribution \f[ p(x)\,{\rm d}x = {\rm e}^{-x}\,{\rm d}x.\f] A simple
        inversion technique is used. */
    double expon();

    /** This function generates a random number from an exponential distribution function with a
        cut-off, defined by the probability distribution \f[ p(x)\,{\rm d}x \propto \begin{cases}
        \;{\rm e}^{-x}\,{\rm d}x &\qquad \text{for $0<x<x_{\text{max}}$,} \\ \;0 &\qquad \text{for
        $x>x_{\text{max}}$.} \end{cases} \f] A simple inversion technique is used. */
    double exponCutoff(double xmax);

    /** This function generates a random direction on the unit sphere, i.e. a couple
        \f$(\theta,\phi)\f$ from the two-dimensional probability density \f[
        p(\theta,\phi)\,d\theta\,d\phi = \left(\frac{\sin\theta}{2}\,d\theta\right)
        \left(\frac{1}{2\pi}\,d\varphi\right).\f] A random direction on the unit sphere can thus be
        constructed by taking two uniform deviates \f${\cal{X}}_1\f$ and \f${\cal{X}}_2\f$, and
        solving the two equations \f[ \begin{split} {\cal{X}}_1 &= \int_0^\theta
        \frac{\sin\theta'\,d\theta'}{2} \\ {\cal{X}}_2 &= \int_0^\varphi \frac{d\varphi'}{2\pi}
        \end{split} \f] for \f$\theta\f$ and \f$\varphi\f$. The solution is readily found, \f[
        \begin{split} \theta &= \arccos\left(2{\cal{X}}_1-1\right) \\ \varphi &= 2\pi\,{\cal {X}}_2.
        \end{split} \f] Once these spherical coordinates are calculated, a Direction object can be
        constructed by calling the constructor Direction::Direction(double theta, double phi). */
    Direction direction();

    /** This function generates a new direction on the unit sphere deviating from a given original
        direction \f$\bf{k}\f$ by a given polar angle \f$\theta\f$ (specified through its cosine)
        and a uniformly random azimuth angle \f$\phi\f$. The function can use an arbitrary
        reference point for the azimuth angle since it is distributed uniformly. We use equation
        (A37) of Bianchi et al. 1996 (ApJ 465,127), which uses the projection of the z-axis as the
        reference point. */
    Direction direction(Direction bfk, double costheta);

    /** This function generates a uniformly distributed random position in a given box (i.e. a
        cuboid lined up with the coordinate axes). */
    Position position(const Box& box);

    //======================== Data Members ========================

private:
    // the state of the random generators; one for each concurrent thread in the simulation
    // (maintaining a separate generator per thread avoids time-consuming data locking)
    vector<vector<unsigned long>> _mtv;
    vector<int> _mtiv;

    // a cached pointer to the ParallelFactory instance associated with this simulation hierarchy
    ParallelFactory* _parfac{nullptr};
};

//////////////////////////////////////////////////////////////////////

#endif

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
class Position;

//////////////////////////////////////////////////////////////////////

/** The Random class offers pseudo-random number generation facilities in a multi-threaded
    environment. At the core of its operation is the function to generate a uniform deviate, i.e. a
    pseudo-random number drawn from the uniform probability distribution over the unit interval.
    Additional functions allow drawing pseudo-random numbers from other frequently-used probability
    distributions.

    To avoid the need for synchronization between multiple execution threads, each thread receives
    its own pseudo-random number generator instance. In fact, for each execution thread in the
    process, the class constructs two thread-local pseudo-random number generators. The \em
    arbitrary generator is seeded with a truly random state obtained from the operating system, so
    that it delivers a unique and unpredictable pseudo-random sequence for each thread, and for
    each execution of the program. The \em predictable generator is always initialized with the
    same fixed state (depending on the value of the user-configurable \em seed property), so that
    it delivers exactly the same pseudo-random sequence for every thread and for every execution of
    the program.

    At any one time, a particular thread in the process uses one of these two generators. Most
    threads employ the arbitrary generator at all times. As an exception, the \em parent thread of
    a Random instance (i.e. the thread that called its setup() function) may sometimes use the
    predictable generator. Specifically, the Random::setupSelfBefore() function (re-)initializes
    the predictable random generator to its fixed initial state (depending on the value of the
    user-configurable \em seed property), and causes the calling thread (i.e., its parent thread)
    to use the predictable generator. The Random::switchToArbitrary() and
    Random::switchToPredictable() functions can be used to switch the parent thread between both
    generators.

    This complicated mechanism is intended to support the following use cases:
      - Repeated execution in single thread/single process mode (with the same \em seed value)
        must produce identical results. This is important, for example, for automated testing.
      - The user can force a different randomization, even in single thread/single process mode,
        by specifying another \em seed value.
      - With parallel threads, tasks are performed in unpredictable order and thus the random
        sequence is interrogated in an unpredictable order; the sequence can just as well be
        unpredictable to begin with.
      - When distributing tasks over multiple processes, each of the threads requires a different
        random sequence.
      - When duplicating tasks in multiple processes, each of the processes must produce an
        identical result (using a single thread in each process) and thus those threads require
        the same random sequence.

    The recommended use of the Random class is to include a single instance in each simulation
    run-time hierarchy. The Parallel subclasses returned from a ParallelFactory in the same
    run-time hierarchy will then call the switchToArbitrary() and switchToPredictable() functions
    at the appropriate times.

    All rng's used in this class are based on the 64-bit Mersenne twister, which offers a
    sufficiently long period and acceptable spectral properties for most purposes. */
class Random : public SimulationItem
{
    ITEM_CONCRETE(Random, SimulationItem, "the default random generator")

    PROPERTY_INT(seed, "the seed for the random generator")
        ATTRIBUTE_MIN_VALUE(seed, "0")
        ATTRIBUTE_MAX_VALUE(seed, "1000000")
        ATTRIBUTE_DEFAULT_VALUE(seed, "0")
        ATTRIBUTE_SILENT(seed)

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function (re-)initializes the predictable random generator to its fixed initial state
        (depending on the value of the user-configurable \em seed property), and causes the calling
        thread (i.e., the parent thread) to use the predictable generator as opposed to the default
        arbitrary generator. The switchToArbitrary() and switchToPredictable() functions can be used
        to switch the parent thread between both generators. */
    void setupSelfBefore() override;

public:
    /** This function switches the calling thread to the arbitrary generator. It should be called
        only from the parent thread, i.e. the thread that called setup() on this instance. */
    void switchToArbitrary();

    /** This function switches the calling thread to the predictable generator. It should be called
        only from the parent thread, i.e. the thread that called setup() on this instance. */
    void switchToPredictable();

    //======================== Other Functions =======================

public:
    /** This function generates a uniform deviate, i.e. a random double precision number in the
        open interval (0,1). The interval borders zero and one are never returned. */
    double uniform();

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

    /** This function generates a random number drawn from an arbitrary probability distribution
        \f$p(x)\,{\text{d}}x\f$ with corresponding cumulative distribution function \f$P(x)\f$. The
        function accepts a discretized version \f$P_i\f$ of the cdf sampled at a set of \f$N\f$
        points \f$x_i\f$. A uniform deviate \f$\cal{X}\f$ is generated, and the equation
        \f${\cal{X}}=P(x)\f$ is solved using linear interpolation on both the axis and probability
        density values. */
    double cdfLinLin(const Array& xv, const Array& Pv);

    /** This function generates a random number drawn from an arbitrary probability distribution
        \f$p(x)\,{\text{d}}x\f$ with corresponding cumulative distribution function \f$P(x)\f$. The
        function accepts a discretized version \f$P_i\f$ of the cdf sampled at a set of \f$N\f$
        points \f$x_i\f$. A uniform deviate \f$\cal{X}\f$ is generated, and the equation
        \f${\cal{X}}=P(x)\f$ is solved using logarithmic interpolation on both the axis and
        probability density values. */
    double cdfLogLog(const Array& xv, const Array& pv, const Array& Pv);
};

//////////////////////////////////////////////////////////////////////

#endif

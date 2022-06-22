/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
class Box;
class Direction;
class Position;
class Vec;

//////////////////////////////////////////////////////////////////////

/** The Random class offers pseudo-random number generation facilities in a multi-threaded
    environment. At the core of its operation is the function to generate a uniform deviate, i.e. a
    pseudo-random number drawn from the uniform probability distribution over the unit interval.
    Additional functions allow drawing pseudo-random numbers from other frequently-used probability
    distributions.

    To avoid the need for synchronization between multiple execution threads, each thread receives
    its own thread-local pseudo-random number generator instance. The parent thread in each process
    (i.e. the thread that calls the setup() function) receives a \em predictable generator. This
    generator is initialized with a fixed state depending only on the value of the
    user-configurable \em seed property, so that it delivers exactly the same pseudo-random
    sequence for every parent thread in every process and for every execution of the program. All
    other (child) threads receive an \em arbitrary generator. These generators are seeded with a
    truly random state obtained from the operating system, so that they deliver a unique and
    unpredictable pseudo-random sequence for each thread, and for each execution of the program.

    This mechanism is intended to support the following use cases. In single thread/single process
    mode, repeated execution (with the same \em seed value) must produce identical results. This is
    important, for example, for automated testing. The user can force a different randomization by
    specifying another \em seed value. In multi-threading and/or multi-processing mode, for tasks
    employing parallelization, each of the threads in each of the processes requires a different
    random sequence. On the other hand, in this mode, for serialized tasks, the (single) thread
    employed in each of the processes must receive the same random sequence. This is important for
    functions that rely on randomness but need to produce the same result in every parallel
    process.

    The recommended use of the Random class is to include a single instance in each simulation
    run-time hierarchy. The use cases discussed above can be achieved as follows. In single
    thread/single process mode, perform all tasks in the parent thread. In multi-threading and/or
    multi-processing mode, perform serial tasks in the parent thread of each process, and perform
    all parallized tasks in a child thread.

    As an additional service, this class allows temporarily installing a predictable random number
    generator with a given seed in the current execution thread through the push() and pop()
    functions. This supports the use case where, usually during setup, the same pseudo-random
    sequence is required in multiple places.

    All random number generators used in this class are based on the 64-bit Mersenne twister, which
    offers a sufficiently long period and acceptable spectral properties for most purposes. */
class Random : public SimulationItem
{
    ITEM_CONCRETE(Random, SimulationItem, "the default random generator")

        PROPERTY_INT(seed, "the seed for the random generator")
        ATTRIBUTE_MIN_VALUE(seed, "0")
        ATTRIBUTE_MAX_VALUE(seed, "1000000")
        ATTRIBUTE_DEFAULT_VALUE(seed, "0")
        ATTRIBUTE_DISPLAYED_IF(seed, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function initializes the random generator for the calling thread (deemed the \em
        parent thread) to its fixed initial state depending on the value of the user-configurable
        \em seed property. */
    void setupSelfBefore() override;

    //=================== Generating random numbers ===================

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

    /** This function generates a random velocity from a three-dimensional Maxwell-Boltzmann
        distribution with velocity dispersion 1, which is equivalent to a Gaussian distribution
        with mean 0 and standard deviation 1 for each of the Cartesian velocity components. The
        implementation simply calls the function gauss() to obtain each of these three velocity
        components. */
    Vec maxwell();

    /** This function generates a random number drawn from an arbitrary probability distribution
        \f$p(x)\,{\text{d}}x\f$ with corresponding cumulative distribution function \f$P(x)\f$. The
        function accepts a discretized version \f$P_i\f$ of the cdf sampled at a set of \f$N\f$
        points \f$x_i\f$. A uniform deviate \f$\cal{X}\f$ is generated, and the equation
        \f${\cal{X}}=P(x)\f$ is solved using linear interpolation (i.e. assuming piece-wise linear
        behavior of the cdf (and equivalently, of the underlying pdf). */
    double cdfLinLin(const Array& xv, const Array& Pv);

    /** This function generates a random number drawn from an arbitrary probability distribution
        \f$p(x)\,{\text{d}}x\f$ with corresponding cumulative distribution function \f$P(x)\f$. The
        function accepts discretized versions \f$p_i\f$ and \f$P_i\f$ of the pdf and cdf sampled at
        a set of \f$N\f$ points \f$x_i\f$. A uniform deviate \f$\cal{X}\f$ is generated, and the
        equation \f${\cal{X}}=P(x)\f$ is solved using a interpolation mechanism that assumes that
        the pdf is linear in log-log space between any two grid points (equivalent to power-law
        behavior), as described below.

        Consider the pdf values \f$p_i\f$ and \f$p_{i+1}\f$ at two consecutive grid points
        \f$x_i\f$ and \f$x_{i+1}\f$. Assuming power-law behavior, the pdf between these two grid
        points can be written as \f[ p(x) = p_i \left(\frac{x}{x_i}\right)^{\alpha_i},
        \quad\mathrm{with}\; \alpha_i = \frac{\ln p_{i+1}/\ln p_i}{\ln x_{i+1}/\ln x_i} \f]

        With \f$\mathcal{X}\f$ a random deviate for which \f$x_i\leq\mathcal{X}\leq x_{i+1}\f$
        happens to be true, we thus need to invert the relation \f[ \mathcal{X}-x_i = \int_{x_i}^x
        p(x')\,\mathrm{d}x' = \int_{x_i}^x p_i \left(\frac{x'}{x_i}\right)^{\alpha_i}\mathrm{d}x' =
        p_i x_i \;\mathrm{gln}\left(-\alpha_i, \frac{x}{x_i}\right) \f] which leads to
        \f[x = x_i \;\mathrm{gexp}\left(-\alpha_i, \frac{\mathcal{X}-x_i}{p_i x_i}\right)\f]
        where \f$\mathrm{gln}(a,x)\f$ is the generalized logarithm and \f$\mathrm{gexp}(a,x)\f$ the
        generalized exponential, defined in the description of respectively the
        SpecialFunctions::gln() and SpecialFunctions::gexp() functions. */
    double cdfLogLog(const Array& xv, const Array& pv, const Array& Pv);

    //=================== Installing a temporary generator ===================

public:
    /** This function pushes the active random number generator for the current thread onto the
        stack for the current thread and establishes a new predictable random number generator,
        initialized with the specified seed, as the active random number generator for the current
        thread. */
    void push(int seed);

    /** This function pops the most recently pushed random number generator from the stack for the
        current thread and establishes it as the active random number generator for the current
        thread. If the stack does not contain a random number generator, the behavior of this
        function is undefined. */
    void pop();
};

//////////////////////////////////////////////////////////////////////

#endif

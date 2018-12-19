/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOCHASTICDUSTEMISSIVITY_HPP
#define STOCHASTICDUSTEMISSIVITY_HPP

#include "DustEmissivity.hpp"
#include <map>
class MultiGrainDustMix;
class SDE_Calculator;
class SDE_Grid;

//////////////////////////////////////////////////////////////////////

/** The StochasticDustEmissivity class calculates the emissivity of a particular dust mix in a
    given radiation field taking into account the temperature distribution of stochastically heated
    dust grains (i.e. without assuming that the dust grains are in local thermal equilibrium) where
    appropriate.

    This calculation requires access to enthalpy data for the various dust grain materials in a
    dust mix, and assumes that the grain size distribution is discretized into a sufficient number
    of size bins. Therefore, when this class is configured, all dust mixes in the media system must
    offer the MultiGrainDustMix interface.

    Two distinct wavelength discretizations (grids) play a role in the calculations performed by
    this class. The radiation field surrounding the dust mix is specified on the simulation's
    radiation field wavelength grid, i.e. the grid returned by the
    Configuration::radiationFieldWLG() function. On the other hand, the resulting dust emissivity
    must be calculated on the simulation's dust emission wavelength grid, i.e. the grid returned by
    the Configuration::dustEmissionWLG() function. In the context of this class, these grids are
    called the input and output wavelength grid, respectively.

    Using the discretization of the dust composition and size distribution into a range of
    representative grains (size bins) and the output wavelength grid described above, in addition
    to a specialized temperature grid constructed by this class, the emissivity in an interstellar
    radiation field \f$J_\lambda\f$ can be written as \f[ \varepsilon_\lambda = \frac{1}{\mu}
    \sum_{b=0}^{N_{\text{size}}-1} \varsigma_{\lambda,b}^{\text{abs}}\,
    \sum_{i=0}^{N_{\text{temp}}-1} P_{b,i}\, B_\lambda(T_i) \f] with \f$\mu\f$ the total dust mass
    of the dust mix, \f$\varsigma_{\lambda,b}^{\text{abs}}\f$ the absorption cross section of the
    \f$b\f$'th size bin, \f$T_i\f$ the \f$i\f$'th temperature grid point, and \f$P_{b,i}\f$ the
    probability of finding a grain of the \f$b\f$'th size bin in the \f$i\f$'th temperature bin.

    The probabilities \f$P_{b,i}\f$ are calculated following a scheme based on Guhathakurta &
    Draine (ApJ 1989), Draine & Li (ApJ 2001), Kruegel (book, 2003), and Misselt et al. (arXiv,
    2008). The scheme is also nicely described by Verstappen (PhD thesis, 2013, section 2.5) and is
    summarized by Camps et al. (A\&A 2015). An overview is presented below. To simplify the
    notation we focus on a single representative dust grain, dropping the index \f$b\f$. Also note
    that the value of the radiation field is now obtained from the input wavelength grid.

    We define a transition matrix \f$A_{f,i}\f$ describing the probability per unit time for a
    grain to transfer from initial temperature bin \f$i\f$ to final temperature bin \f$f\f$. The
    transition matrix elements in the case of heating \f$(f>i)\f$ are given by \f[ A_{f,i} = 4\pi\,
    \varsigma_{\lambda_{fi}}^{\text{abs}}\,J_{\lambda_{fi}}\, \frac{hc\,\Delta H_f}{(H_f-H_i)^3}
    \f] where \f$H_f\f$ and \f$H_i\f$ are the enthalpies of the final and initial temperature bins,
    \f$\Delta H_f\f$ is the width of the final temperature bin, and \f$\lambda_{fi}\f$ is the
    transition wavelength which can be obtained from \f[ \lambda_{fi}=\frac{hc}{H_f-H_i}. \f] We
    assume that cooling transitions occur only to the next lower level, so that \f$A_{f,i}=0\f$ for
    \f$f<i-1\f$ and \f[ A_{i-1,i} = \frac{4\pi}{H_i-H_{i-1}}\, \int_0^\infty
    \varsigma_{\lambda}^{\text{abs}}\, B_{\lambda}(T_i)\, {\text{d}}\lambda. \f] The diagonal
    matrix elements are defined as \f[ A_{i,i} = -\sum_{f\ne i} A_{f,i}\f] however as we will see
    below there is no need to explicitly calculate these values.

    Assuming a steady state situation, the probabilities \f$P_{i}\f$ can be obtained from the
    transition matrix by solving the set of \f$N\f$ linear equations \f[ \sum_{i=0}^{N-1} A_{f,i}
    \,P_i=0 \qquad f=0,...,N-1 \f] along with the normalization condition \f[ \sum_{i=0}^{N-1} P_i
    = 1, \f] where \f$N\f$ is the number of temperature bins. Because the matrix values for
    \f$f<i-1\f$ are zero these equations can be solved by a recursive procedure of computational
    order \f${\mathcal{O}(N^2)}\f$. To avoid numerical instabilities caused by the negative
    diagonal elements, the procedure employs a well-chosen linear combination of the original
    equations. This leads to the following recursion relations for the adjusted matrix elements
    \f$B_{f,i}\f$, the unnormalized probabilities \f$X_i\f$, and finally the normalized
    probabilities \f$P_{i}\f$: \f{align*} B_{N-1,i} &= A_{N-1,i} & i=0,\ldots,N-2 \\ B_{f,i} &=
    B_{f+1,i}+A_{f,i} & f=N-2,\ldots,1;\,i=0,\ldots,f-1 \\ X_0 &= 1 \\ X_i &=
    \frac{\sum_{j=0}^{i-1}B_{i,j}X_j}{A_{i-1,i}} & i=1,\ldots,N-1 \\ P_i &=
    \frac{X_i}{\sum_{j=0}^{N-1}X_j} & i=0,\ldots,N-1 \f}

    */
class StochasticDustEmissivity : public DustEmissivity
{
    ITEM_CONCRETE(StochasticDustEmissivity, DustEmissivity, "Stochastic (Non-LTE) dust emissivity calculator")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor destructs the data structures created during setup. */
    ~StochasticDustEmissivity();

    /** This function precalculates and stores a bunch of data for each representative dust grain
        (size bin) in the simulation with the intention to accelerate the actual dust emissivity
        computations. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the dust emissivity spectrum \f$\varepsilon_{\ell'}\f$ for the
        specified dust mix residing in a radiation field with the specified mean intensities
        \f$J_\ell\f$. The input and output arrays are discretized on the wavelength grids returned
        by the Configuration::radiationFieldWLG() and Configuration::dustEmissionWLG() functions,
        repectively.

        If the specified material mix does not offer the MultiGrainDustMix interface, the function
        returns zero emissivity. */
    Array emissivity(const MaterialMix* mix, const Array& Jv) const override;

    //========================= Data members =======================

private:
    // the size of the output wavelength grid
    int _numOutWavelengths{0};

    // setupSelfBefore adds all temperature grids to this list so that they stay around and can be destructed
    vector<const SDE_Grid*> _grids;

    // setupSelfBefore adds a calculator to these maps for each representative dust grain (size bin) in the simulation
    std::map<std::pair<const MultiGrainDustMix*,int>, const SDE_Calculator*> _calculatorsA;     // coarse grid
    std::map<std::pair<const MultiGrainDustMix*,int>, const SDE_Calculator*> _calculatorsB;     // medium grid
    std::map<std::pair<const MultiGrainDustMix*,int>, const SDE_Calculator*> _calculatorsC;     // fine grid
};

////////////////////////////////////////////////////////////////////

#endif

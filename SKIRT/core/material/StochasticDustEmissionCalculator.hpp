/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOCHASTICDUSTEMISSIONCALCULATOR_HPP
#define STOCHASTICDUSTEMISSIONCALCULATOR_HPP

#include "Array.hpp"
#include "StoredTable.hpp"
class SimulationItem;
class SDE_Calculator;
class SDE_TemperatureGrid;

////////////////////////////////////////////////////////////////////

/** StochasticDustEmissionCalculator is a helper class to calculate the combined emissivity
    spectrum \f$\varepsilon_{\lambda}\f$ of a set of representative dust grain populations embedded
    in a given radiation field \f$J_\lambda\f$, taking into account the temperature distribution of
    stochastically heated dust grains (i.e. without assuming that the dust grains are in local
    thermal equilibrium) where appropriate.

    The class handles multiple representative grains that together describe a dust mixture. These
    are called \em bins because the calculator is typically used to handle the various grain size
    bins in a dust mix. The calculator assumes that, for each material type, the size distribution
    is discretized into a sufficient number of size bins. The calculator also requires access to
    enthalpy data for the various dust grain materials in the dust mix.

    A client of the class must first call the precalculate() function for each bin to supply the
    absorption cross sections and enthalpy data for the representative grain population
    corresponding to that bin. The emissivity() function can then be used to obtain the combined
    emissivity spectrum for the representative grains in all bins. The embedding radiation field is
    specified by the mean intensities \f$(J_\lambda)_k\f$, which must be discretized on the
    simulation's radiation field wavelength grid as returned by the
    Configuration::radiationFieldWLG() function. The emissivity spectrum
    \f$(\varepsilon_\lambda)_\ell\f$ produced by the second function is discretized on the
    wavelength grid returned by the Configuration::dustEmissionWLG() function.

    Using the discretization of the dust composition and size distribution into a range of
    representative grains (size bins) and the output wavelength grid described above, in addition
    to a specialized temperature grid constructed by this class, the emissivity in an interstellar
    radiation field \f$J_\lambda\f$ can be written as \f[ \varepsilon_\lambda =
    \sum_{b=0}^{N_{\text{size}}-1} \varsigma_{\lambda,b}^{\text{abs}}\,
    \sum_{i=0}^{N_{\text{temp}}-1} P_{b,i}\, B_\lambda(T_i) \f] with
    \f$\varsigma_{\lambda,b}^{\text{abs}}\f$ the absorption cross section of the \f$b\f$'th size
    bin, \f$T_i\f$ the \f$i\f$'th temperature grid point, and \f$P_{b,i}\f$ the probability of
    finding a grain of the \f$b\f$'th size bin in the \f$i\f$'th temperature bin.

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
    transition wavelength which can be obtained from \f[ \lambda_{fi}=\frac{hc}{H_f-H_i}. \f] If
    the simulation's configuration requests the inclusion of the cosmic microwave background (CMB)
    as an additional source for dust heating, the corresponding term is added to the heating matrix
    element, similar to the treatment of the energy balance equation in the
    EquilibriumDustEmissionCalculator class. We assume that cooling transitions occur only to the
    next lower level, so that \f$A_{f,i}=0\f$ for \f$f<i-1\f$ and \f[ A_{i-1,i} =
    \frac{4\pi}{H_i-H_{i-1}}\, \int_0^\infty \varsigma_{\lambda}^{\text{abs}}\, B_{\lambda}(T_i)\,
    {\text{d}}\lambda. \f] The diagonal matrix elements are defined as \f[ A_{i,i} = -\sum_{f\ne i}
    A_{f,i}\f] however as we will see below there is no need to explicitly calculate these values.

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
class StochasticDustEmissionCalculator
{
public:
    /** The destructor destructs the data structures allocated by the precalculate() function. */
    ~StochasticDustEmissionCalculator();

    /** This function precalculates and stores information used to calculate the emissivity
        spectrum for a particular bin (i.e. representative grain) to be handled by the calculator.
        It must be called once for each bin.

        When it is first called, the function uses its first argument to obtain the simulation's
        radiation field and dust emission wavelength grids, and it builds a set of temperature
        grids for use in the calculator: a coarse grid for quickly determining the appropriate
        temperature range, and medium and fine grids for performing the actual probability
        calculations.

        The second and third function arguments specify the absorption cross sections
        \f$\varsigma^\text{abs}_i\f$ for the representative grain corresponding to the current bin
        on some fine wavelength grid \f$\lambda_i\f$. The remaining arguments specify an identifier
        for the type of grain material in the bin, the bulk density of the grain material, the
        average mass of a dust grain in the current bin, and a stored table containing the specific
        enthalpy (per unit volume) for the grain material as a function of temperature.

        The function stores the absorption cross sections interpolated on the radiation field and
        dust emission wavelength grids and it precalculates Planck-integrated absorption cross
        sections on each of the constructed temperature grids through integration over the fine
        wavelength grid specified as an argument. Furthermore, the function precalculates the
        heating and cooling rates used for the stochastic probability calculations, barring the
        input radiation field dependency, again on each of the constructed temperature grids. */
    void precalculate(SimulationItem* item, const Array& lambdav, const Array& sigmaabsv, string grainType,
                      double bulkDensity, double meanMass, const StoredTable<1>& enthalpy);

    /** This function returns the size of the memory, in bytes, allocated by the precalculate()
        function so far. This information can be used for logging purposes. */
    size_t allocatedBytes() const;

    /** This function returns the emissivity spectrum per hydrogen atom
        \f$(\varepsilon_\lambda)_\ell\f$ of the dust mix (or rather of the corresponding mixture of
        representative grain populations) when embedded in the radiation field specified by the
        mean intensities \f$(J_\lambda)_k\f$, taking into account the temperature distribution of
        stochastically heated dust grains. The input and output arrays are discretized on the
        wavelength grids returned by the Configuration::radiationFieldWLG() and
        Configuration::dustEmissionWLG() functions, repectively. If the precalculate() function has
        not been called for at least one bin, the behavior of this function is undefined. */
    Array emissivity(const Array& Jv) const;

    //======================== Data Members ========================

private:
    // wavelength grids and CMB
    Array _rflambdav;   // radiation field wavelength grid (RFWLG) -- indexed on k
    Array _rfdlambdav;  // radiation field wavelength grid bin widths -- indexed on k
    Array _Bcmbv;       // cosmic microwave background radiation field, if needed -- indexed on k
    Array _emlambdav;   // dust emission wavelength grid (EMWLG) -- indexed on ell

    // temperature grids
    const SDE_TemperatureGrid* _gridA{nullptr};  // coarse grid
    const SDE_TemperatureGrid* _gridB{nullptr};  // medium grid
    const SDE_TemperatureGrid* _gridC{nullptr};  // fine grid

    // calculators for each representative dust grain (size bin) -- indexed on b
    vector<const SDE_Calculator*> _calculatorsA;  // coarse grid
    vector<const SDE_Calculator*> _calculatorsB;  // medium grid
    vector<const SDE_Calculator*> _calculatorsC;  // fine grid

    // other properties for each representative dust grain (size bin) -- indexed on b
    vector<string> _grainTypes;        // the grain type identifier
    vector<double> _meanMasses;        // mean mass of a grain
    vector<double> _maxEnthalpyTemps;  // maximum temperature for the enthalpy data
};

////////////////////////////////////////////////////////////////////

#endif

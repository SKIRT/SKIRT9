/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EQUILIBRIUMDUSTEMISSIONCALCULATOR_HPP
#define EQUILIBRIUMDUSTEMISSIONCALCULATOR_HPP

#include "Array.hpp"
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** EquilibriumDustEmissionCalculator is a helper class to calculate the equilibrium temperature
    \f$T_{\text{eq}}\f$ and the emissivity spectrum \f$\varepsilon_{\lambda}\f$ of representative
    dust grain populations embedded in a given radiation field \f$J_\lambda\f$.

    The class is equipped to handle multiple (in principle independent) representative grains.
    These are called \em bins because the calculator is typically used to handle the various grain
    size bins in a dust mix. Indeed, the equilibrium temperature and the emissivity spectrum are
    nonlinear functions of the grain properties, and thus a single representative grain cannot
    usually accurately represent a dust mix.

    A client of the class must first call the precalculate() function for each bin (in order of bin
    index \f$b\f$) to supply the absorption cross sections for the representative grain population
    corresponding to that bin. The equilibriumTemperature() function can then be used to obtain the
    equilibrium temperature for a given bin, and the emissivity() function can be used to obtain
    the combined emissivity spectrum for the representative grains in all bins. In both cases, the
    embedding radiation field is specified by the mean intensities \f$(J_\lambda)_k\f$, which must
    be discretized on the simulation's radiation field wavelength grid as returned by the
    Configuration::radiationFieldWLG() function. The emissivity spectrum
    \f$(\varepsilon_\lambda)_\ell\f$ produced by the second function is discretized on the
    wavelength grid returned by the Configuration::dustEmissionWLG() function.

    The equilibrium temperature \f$T_{\text{eq},b}\f$ for bin with index \f$b\f$ is obtained from
    the energy balance equation, \f[ \int_0^\infty \varsigma_{\lambda,b}^{\text{abs}}\, J_\lambda\,
    {\text{d}}\lambda = \int_0^\infty \varsigma_{\lambda,b}^{\text{abs}}\,
    B_\lambda(T_{\text{eq},b})\, {\text{d}}\lambda \f] where \f$B_\lambda(T)\f$ is the Planck
    function. The left-hand side is integrated over the radiation field wavelength grid, and the
    right-hand side is precalculated for a range of temperatures through integration over the
    wavelength grid and cross sections passed to the precalculate() function. If the simulation's
    configuration requests the inclusion of the cosmic microwave background (CMB) as an additional
    source for dust heating, the corresponding term is added to the left-hand side of the energy
    balance equation: \f[ \int_0^\infty \varsigma_{\lambda,b}^{\text{abs}}\, \left[
    B_\lambda\left((1+z) T_\mathrm{CMB}^{z=0}\right) + J_\lambda\right] \, {\text{d}}\lambda =
    \int_0^\infty \varsigma_{\lambda,b}^{\text{abs}}\, B_\lambda(T_{\text{eq},b})\,
    {\text{d}}\lambda \f] where \f$z\f$ is the redshift at which the simulated model resides and
    \f$T_\mathrm{CMB}^{z=0} = 2.725\,\mathrm{K}\f$.

    The equilibrium emissivity spectrum of all bins combined embedded in a radiation field
    \f$J_\lambda\f$ can then be written as \f[ \varepsilon_\lambda = \sum_{b=0}^{N_{\text{bins}}-1}
    \varsigma_{\lambda,b}^{\text{abs}}\, B_\lambda(T_{\text{eq},b}) \f] with
    \f$\varsigma_{\lambda,b}^{\text{abs}}\f$ the absorption cross section of the \f$b\f$'th
    representative grain and \f$T_{\text{eq},b}\f$ the equilibrium temperature of that grain. */
class EquilibriumDustEmissionCalculator
{
public:
    /** This function precalculates and stores information used to calculate the equilibrium
        temperature and emissivity spectrum for a particular bin (i.e. representative grain) to be
        handled by the calculator. It must be called once for each bin; the order of the calls
        establishes the order of the bin index \f$b\f$.

        When it is first called, the function uses its first argument to obtain the simulation's
        radiation field and dust emission wavelength grids, and it builds a temperature grid for
        use in the calculator.

        The second and third function arguments specify the absorption cross sections
        \f$\varsigma^\text{abs}_i\f$ for the representative grain corresponding to the current bin
        on some fine wavelength grid \f$\lambda_i\f$. The function stores the absorption cross
        sections interpolated on the radiation field and dust emission wavelength grids and it
        precalculates Planck-integrated absorption cross sections on an appropriate temperature
        grid through integration over the fine wavelength grid specified as an argument. */
    void precalculate(SimulationItem* item, const Array& lambdav, const Array& sigmaabsv);

    /** This function returns the size of the memory, in bytes, allocated by the precalculate()
        function so far. This information can be used for logging purposes. */
    size_t allocatedBytes() const;

    /** This function returns the number of bins (i.e. representative grains) added by the
        precalculate() function so far. */
    int numBins() const;

    /** This function returns the equilibrium temperature \f$T_{\text{eq},b}\f$ of the
        representative grain corresponding to the bin with specified index \f$b\f$ when embedded in
        the radiation field specified by the mean intensities \f$(J_\lambda)_k\f$, which must be
        discretized on the simulation's radiation field wavelength grid as returned by the
        Configuration::radiationFieldWLG() function. If the precalculate() function has not been
        called for the specified bin, the behavior of this function is undefined. */
    double equilibriumTemperature(int b, const Array& Jv) const;

    /** This function returns the emissivity spectrum per hydrogen atom
        \f$(\varepsilon_\lambda)_\ell\f$ of the dust mix (or rather of the corresponding mixture of
        representative grain populations) when embedded in the radiation field specified by the
        mean intensities \f$(J_\lambda)_k\f$, assuming that the dust grains are in local thermal
        equilibrium. The input and output arrays are discretized on the wavelength grids returned
        by the Configuration::radiationFieldWLG() and Configuration::dustEmissionWLG() functions,
        repectively. If the precalculate() function has not been called for at least one bin, the
        behavior of this function is undefined. */
    Array emissivity(const Array& Jv) const;

    //======================== Data Members ========================

private:
    Array _rflambdav;   // radiation field wavelength grid (RFWLG) -- indexed on k
    Array _rfdlambdav;  // radiation field wavelength grid bin widths -- indexed on k
    Array _Bcmbv;       // cosmic microwave background radiation field, or zeroes -- indexed on k
    Array _emlambdav;   // dust emission wavelength grid (EMWLG) -- indexed on ell
    Array _Tv;          // temperature grid for the integrated absorption cross sections -- indexed on p

    vector<Array> _rfsigmaabsvv;  // absorption cross sections on the RFWLG for each bin -- indexed on b,k
    vector<Array> _emsigmaabsvv;  // absorption cross sections on the EMWLG for each bin -- indexed on b,ell
    vector<Array> _planckabsvv;   // Planck-integrated absorption cross sections for each bin -- indexed on b,p
};

////////////////////////////////////////////////////////////////////

#endif

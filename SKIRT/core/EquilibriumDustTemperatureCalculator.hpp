/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef EQUILIBRIUMDUSTTEMPERATURECALCULATOR_HPP
#define EQUILIBRIUMDUSTTEMPERATURECALCULATOR_HPP

#include "Array.hpp"
class SimulationItem;
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** EquilibriumDustTemperatureCalculator is a helper class to calculate the equilibrium temperature
    \f$T_{\text{eq}}\f$ of a representative dust grain population embedded in a given radiation
    field. The class is equipped to handle multiple (in principle independent) representative
    grains. These are called \em bins because the calculator is typically used to handle the
    various grain size bins in a dust mix. Indeed, the equilibrium temperature is a nonlinear
    function of the grain properties, and thus a single representative grain cannot usually
    accurately represent a dust mix. Refer to the description of the MultiGrainDustMix class for
    more information.

    A client of the class must first call the precalculate() function for each bin (in order of bin
    index \f$b\f$) to supply the absorption cross sections for the representative grain population
    corresponding to that bin. The equilibriumTemperature() function can then be used to obtain the
    equilibrium temperature for a given bin. The embedding radiation field is specified by the mean
    intensities \f$(J_\lambda)_\ell\f$, which must be discretized on the simulation's radiation
    field wavelength grid as returned by the Configuration::radiationFieldWLG() function.

    The equilibrium temperature is obtained from the energy balance equation, \f[ \int_0^\infty
    \varsigma^\text{abs}(\lambda) \,J_\lambda(\lambda) \,\text{d}\lambda = \int_0^\infty
    \varsigma^\text{abs}(\lambda) \,B_\lambda(T_\text{eq},\lambda) \,\text{d}\lambda, \f] where the
    left-hand side is integrated over the radiation field wavelength grid, and the right-hand side
    is precalculated for a range of temperatures through integration over the wavelength grid and
    cross sections passed to the precalculate() function. */
class EquilibriumDustTemperatureCalculator
{
public:
    /** This function precalculates and stores information used to calculate the equilibrium
        temperature for a particular bin (i.e. representative grain) to be handled by the
        calculator. It must be called once for each bin; the order of the calls establishes the
        order of the bin index \f$b\f$.

        When it is first called, the function obtains (a pointer to) the simulation's radiation
        field wavelength grid, and builds a temperature grid for use in the calculator. The first
        argument is used to retrieve the radiation field wavelength grid from the simulation's
        configuration.

        The second and third arguments specify the absorption cross sections
        \f$\varsigma^\text{abs}_\ell\f$ for the representative grain corresponding to the current
        bin on some fine wavelength grid \f$\lambda_\ell\f$. The function stores the absorption
        cross sections interpolated on the radiation field wavelength grid (to facilitate the
        calculation of the input side of the energy balance equation) and it precalculates
        Planck-integrated absorption cross sections on an appropriate temperature grid
        (corresponding to the output side of the energy balance equation) through integration over
        the fine wavelength grid specified as an argument. */
    void precalculate(SimulationItem* item, const Array& lambdav, const Array& sigmaabsv);

    /** This function returns the size of the memory, in bytes, allocated by the precalculate()
        function so far. This information can be used for logging puroposes. */
    size_t allocatedBytes() const;

    /** This function returns the number of bins (i.e. representative grains) added by the
        precalculate() function so far. */
    int numBins() const;

    /** This function returns the equilibrium temperature \f$T_{\text{eq}}\f$ of the representative
        grain corresponding to the bin with specified index \f$b\f$ when embedded in the radiation
        field specified by the mean intensities \f$(J_\lambda)_\ell\f$, which must be discretized
        on the simulation's radiation field wavelength grid as returned by the
        Configuration::radiationFieldWLG() function. If the precalculate() function has not been
        called for the specified bin, the behavior of this function is undefined. */
    double equilibriumTemperature(int b, const Array& Jv) const;

    //======================== Data Members ========================

private:
    WavelengthGrid* _radiationFieldWLG{nullptr};  // the radiation field wavelength grid (RFWLG) -- indexed on ell
    Array _Tv;                   // temperature grid for the integrated absorption cross sections -- indexed on p
    vector<Array> _rfsigmaabsvv; // absorption cross sections on the RFWLG for each bin -- indexed on b,ell
    vector<Array> _planckabsvv;  // Planck-integrated absorption cross sections for each bin -- indexed on b,p
};

////////////////////////////////////////////////////////////////////

#endif

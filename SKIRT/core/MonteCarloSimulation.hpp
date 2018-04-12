/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP

#include "Simulation.hpp"

//////////////////////////////////////////////////////////////////////

/** The MonteCarloSimulation class is the general abstract base class describing Monte Carlo
    simulations. Running a Monte Carlo simulation with SKIRT essentially comes down to constructing
    an instance of one of the subclasses of the MonteCarloSimulation base class and invoking the
    setupAndRun() function on it. The general MonteCarloSimulation class manages ... . */
class MonteCarloSimulation : public Simulation
{
    ITEM_CONCRETE(MonteCarloSimulation, Simulation, "a Monte Carlo simulation")
        ATTRIBUTE_SUB_PROPERTIES_FIRST(MonteCarloSimulation)

    PROPERTY_DOUBLE(numPackages, "the total number of photon packages")
        ATTRIBUTE_MIN_VALUE(numPackages, "[0")
        ATTRIBUTE_MAX_VALUE(numPackages, "1e15]")
        ATTRIBUTE_DEFAULT_VALUE(numPackages, "1e6")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches some frequently used pointers. */
    void setupSelfBefore() override;

    //======== Setters & Getters for Discoverable Attributes =======

    /** \fn numPackages
        Photon packages are launched in chunks of the same size. The chunk size is determined
        automatically during setup. The number of photon packages actually launched is always an
        integer multiple of the chunk size, and thus may be slightly above the specified number.
        Unless the specified number of photon packages is exactly equal to zero, a simulation
        always launches at least one chunk.

        The maximum number of photon packages per wavelength is somewhat arbitrarily set to 1e15.
        This function throws an error if a larger number is specified. The argument is of type
        double (which can exactly represent all integers up to 9e15) rather than 64-bit integer to
        avoid implementing yet another discoverable property type. As a side benefit, one can use
        exponential notation to specify a large number of photon packages. */

public:
    /** This function puts the simulation in emulation mode. Specifically, it sets an internal flag
        that can be queried other classes and it sets the number of photon packages to zero. */
    void setEmulationMode();

    /** This function returns true if the simulation has been put in emulation mode. */
    bool emulationMode();

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the simulation, which depends on the (lack of)
        symmetry in the media geometries. A value of 1 means spherical symmetry, 2 means axial
        symmetry and 3 means none of these symmetries. The media component with the least symmetry
        (i.e. the highest dimension) determines the result for the whole simulation. */
    int dimension() const;

protected:
    /** This function actually runs the simulation. */
    void runSelf() override;

    //======================== Data Members ========================

private:
    // *** data member to remember whether emulation mode is enabled
    bool _emulationMode{false};
};

////////////////////////////////////////////////////////////////////

#endif

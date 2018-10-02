/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATIONMODE_HPP
#define SIMULATIONMODE_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** SimulationMode is the abstract base class for a simple class hierarchy that allows the user to
    set the overall operation mode of a SKIRT simulation and configure options for it. Each
    subclass corresponds to a particular fundamental mode of operation, and it offers configuration
    options relevant to that mode. For example, a simulation mode calculating the dust temperature
    self-consistently (taking into account dust self-absorption through iteration) would require
    more extensive configuration options than a mode restricted to calculating the extinction in
    the optical wavelength range.

    In practice, the SimulationMode subclasses bundle configuration options related to the Monte
    Carlo photon lifecycle (e.g. when a photon packet should be terminated), the media state
    (e.g. the wavelength resolution of the radiation field grid), and the iterative simulation
    process (e.g. the convergence criteria).

    The simulation mode should be included as one of the first properties of the top-level
    MonteCarloSimulation class, requiring the user to select a fundamental mode early in the
    configuration process. This in turn facilitates limiting the options offered during the
    remaining configuration process to those relevant for the selected mode. */
class SimulationMode : public SimulationItem
{
    ITEM_ABSTRACT(SimulationMode, SimulationItem, "a simulation mode")
    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif

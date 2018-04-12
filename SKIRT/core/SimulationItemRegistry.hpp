/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATIONITEMREGISTRY_HPP
#define SIMULATIONITEMREGISTRY_HPP

#include "Basics.hpp"
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** The SimulationItemRegistry class manages the registration of all discoverable SimulationItem
    subclasses in the \c SKIRT program (defined in the \c core project subdirectory). A single
    instance of the SimulationItemRegistry class must be constructed just after program startup,
    and certainly before any parallel threads are started. A good place is early in the main()
    function. The program should not use the exit() or abort() functions, but rather let the main()
    function run to normal completion and return an exit code. */
class SimulationItemRegistry final
{
public:
    /** The constructor registers all discoverable SimulationItem subclasses in the \c SKIRT
        program with the item registry, including them in a single SMILE schema called 'SKIRT'. The
        first argument \em version specifies the version of this schema definition. The second
        argument \em format specifies the version of the described data format (listed on the root
        element). The constructor is \em not thread-safe and may be called only during program
        startup from a single thread. */
    SimulationItemRegistry(string version, string format);

    /** This static function returns a pointer to the 'SKIRT' schema definition. Ownership remains
        with the registry. This function is thread-safe and may called at any time after
        construction and before destruction of the SimulationItemRegistry instance. */
    static const SchemaDef* getSchemaDef();

    /** The destructor releases the global memory managed by the item registry. */
    ~SimulationItemRegistry();
};

////////////////////////////////////////////////////////////////////

#endif

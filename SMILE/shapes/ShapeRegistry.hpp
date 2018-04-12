/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHAPEREGISTRY_HPP
#define SHAPEREGISTRY_HPP

#include "Basics.hpp"
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** The ShapeRegistry class manages the registration of all Item subclasses defined and used in the
    'shapes' program. A single instance of the ShapeRegistry class must be constructed just after
    program startup, and certainly before any parallel threads are started. A good place is early
    in the main() function. The program should not use the exit() or abort() functions, but rather
    let the main() function run to normal completion and return an exit code. */
class ShapeRegistry final
{
public:
    /** The constructor registers all Item subclasses defined and used in the 'shapes' program with
        the item registry, including them in a single SMILE schema called 'Shapes'. The function is \em not
        thread-safe and may be called only during program startup from a single thread. */
    ShapeRegistry();

    /** This static function returns a pointer to the 'Shapes' schema definition. Ownership remains
        with the registry. This function is thread-safe and may called at any time after
        construction and before destruction of the ShapeRegistry instance. */
    static const SchemaDef* getSchemaDef();

    /** The destructor releases the global memory managed by the item registry. */
    ~ShapeRegistry();
};

////////////////////////////////////////////////////////////////////

#endif

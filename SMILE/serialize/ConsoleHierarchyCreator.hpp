/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONSOLEHIERARCHYCREATOR_HPP
#define CONSOLEHIERARCHYCREATOR_HPP

#include "Basics.hpp"
class DoubleListPropertyHandler;
class DoublePropertyHandler;
class Item;
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** This class offers a static function to create a SMILE dataset representation in memory through
    user interaction via the console. */
class ConsoleHierarchyCreator final
{
public:
    /** Creates a fresh SMILE dataset for the specified schema definition by asking questions
        (and receiving answers) via the console, and returns a pointer to the root item of the
        hierarchy (handing over ownership for the complete hierarchy to the caller). If the
        hierarchy can't be created due to some error condition the function throws a fatal error.
        */
    static std::unique_ptr<Item> create(const SchemaDef* schema);

    /** Prompts the console user for a double value, using the given message prefix and the
        information provided by the specified property handler. */
    static double promptForDouble(string prefix, const DoublePropertyHandler* handler);

    /** Prompts the console user for a comma-separated list of double values, using the given
        message prefix and the information provided by the specified property handler. */
    static vector<double> promptForDoubleList(string prefix, const DoubleListPropertyHandler* handler);
};

////////////////////////////////////////////////////////////////////

#endif

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XMLHIERARCHYCREATOR_HPP
#define XMLHIERARCHYCREATOR_HPP

#include "Basics.hpp"
class Item;
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** This class offers static functions to load a SMILE dataset from a file or a string containing
    an XML serialization into a memory-based dataset representation. */
class XmlHierarchyCreator final
{
public:
    /** Creates a fresh memory representation of a SMILE dataset for the specified schema
        definition reflecting the XML serialization in the specified file, and returns a pointer to
        the root item of the hierarchy (handing over ownership for the complete hierarchy to the
        caller). If the hierarchy can't be created due to some error condition the function throws
        a fatal error. */
    static std::unique_ptr<Item> readFile(const SchemaDef* schema, string filepath);

    /** Creates a fresh memory representation of a SMILE dataset for the specified schema
        definition reflecting the XML serialization in the specified \em contents string, and
        returns a pointer to the root item of the hierarchy (handing over ownership for the
        complete hierarchy to the caller). The \em description argument provides a human readable
        string to identify the contents string in error messages. If the hierarchy can't be created
        due to some error condition the function throws a fatal error. */
    static std::unique_ptr<Item> readString(const SchemaDef* schema, string contents, string description);
};

////////////////////////////////////////////////////////////////////

#endif

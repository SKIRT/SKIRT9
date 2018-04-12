/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef XMLHIERARCHYWRITER_HPP
#define XMLHIERARCHYWRITER_HPP

#include "Basics.hpp"
class Item;
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** This class offers a static function to write the structure and properties of a SMILE dataset
    representation in memory to an XML file. The XML file contains sufficient information to
    reconstruct a fresh copy of the SMILE dataset. */
class XmlHierarchyWriter final
{
public:
    /** Writes the structure and properties of the specified SMILE dataset described by the given
        schema definition to an XML file with the specified file path. The optional last argument
        specifies a producer identification string to be included as an attribute on the root
        element. If an error occurs, this function throws a fatal error. */
    static void write(Item* item, const SchemaDef* schema, string filePath, string producer = string());
};

////////////////////////////////////////////////////////////////////

#endif

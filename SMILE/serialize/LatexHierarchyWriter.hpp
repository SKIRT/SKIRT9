/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LATEXHIERARCHYWRITER_HPP
#define LATEXHIERARCHYWRITER_HPP

#include "Basics.hpp"
class Item;
class SchemaDef;

////////////////////////////////////////////////////////////////////

/** This class offers a static function to write the structure and properties of a SMILE dataset
    representation in memory to a LaTeX source file so that it can be easily typeset. The class
    supports a number of text and unit string replacements to improve the quality of the TeX
    output. This capability is currently designed rather specifically for the contents of SKIRT
    parameter files. */
class LatexHierarchyWriter final
{
public:
    /** Writes the structure and properties of the specified SMILE dataset described by the given
        schema definition to a LaTeX source file with the specified file path. The \em dataset
        argument describes the dataset being written; e.g. it could be the filename from which the
        dataset has been loaded. The optional last argument specifies a producer identification
        string to be included as an attribute on the root element. If an error occurs, this
        function throws a fatal error. */
    static void write(Item* item, const SchemaDef* schema, string filePath, string dataset, string producer = string());
};

////////////////////////////////////////////////////////////////////

#endif

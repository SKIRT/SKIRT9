/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SHAPESCOMMANDLINEHANDLER_HPP
#define SHAPESCOMMANDLINEHANDLER_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/**
This class offers a static function to handle the command line arguments for the SMILE example
command-line utility 'shapes'. The function parses the command line arguments and calls the
appropriate high-level functions to perform the actions requested by the user. When invoked
with invalid command line arguments, it prints a brief help message. The following command line
arguments are currently supported; their order is not significant.

\verbatim
shapes <shapes_dataset_filepath> [-o <tiff_output_filepath>]
shapes -x <library_dirpath>
\endverbatim

The first form generates a TIFF image based on the instructions provided in the shapes
parameter file (a SMILE dataset) specified as the first command line argument. The \c .shapes
filename extension may be included or omitted. If the optional output filepath is present (the
-o option), it specifies the path to the image file that will be created according to the
instructions in the parameter file. The \c .tiff filename extension is added unless the
filepath already has the \c .tiff filename extension. If the -o option is not present,
the output filepath is derived from the parameter filepath by replacing the \c .smile filename
extension by the .tiff filename extension.

The second form using the -x option causes a SMILE schema file to be generated corresponding
to the metadata in the source code (i.e. in the class definitions of the Item subclasses).
The file is named 'shapes.smile'  and is placed in the specified directory.
*/
class ShapesCommandLineHandler final
{
public:
    /** This function processes the command line arguments and invokes the appropriate high-level
        functions to perform the actions requested by the user. The function returns an appropriate
        program exit value. */
    static int perform();
};

////////////////////////////////////////////////////////////////////

#endif

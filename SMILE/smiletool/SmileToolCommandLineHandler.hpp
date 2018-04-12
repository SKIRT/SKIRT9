/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SMILETOOLCOMMANDLINEHANDLER_HPP
#define SMILETOOLCOMMANDLINEHANDLER_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/**
This class offers a static function to processes the command line arguments for the smiletool
utility and calls the appropriate high-level functions to perform the actions requested by the
user. When invoked with invalid command line arguments, it prints a brief help message. The
following command line arguments are supported; they are all optional and their order is not
significant.

\verbatim
    smiletool [-i <input_dataset_filepath>]
              [-o <output_dataset_filepath>] [-t <output_latex_filepath>]
              [-l <library_dirpath>] [-s <schema_filename>]
\endverbatim

If an input dataset filepath is present (the -i option), the specified dataset is loaded into
memory. Otherwise, an interactive console query and answer session is conducted to create a new
dataset in memory.

If an output dataset filepath is present (the -o option), the memory dataset is saved to the
specified file. Similarly, if an output latex filepath is present (the -t option), the memory
dataset is converted to LaTeX source form and saved to the specified file. In both cases,
an appropriate filename extension is added if needed, and an existing file with the final
name gets overwritten without warning. The -o and -t options can both be specified at the same
time. If neither is specified, the user is prompted for a dataset filename at the console.
In this case, an appropriate filename extension is added if needed as well, however the name
of an existing file will not be accepted.

To operate the tool in batch, i.e. without console prompting, both an input dataset filepath
(the -i option) and at least one output filepath (-o and/or -t options) must be present.

The library directory path (the -l option) specifies the location of SMILE schema files. If
this option is missing, the current directory is used instead. The schema file name (the -s
option) specifies the schema to be used for interpreting or constructing a dataset. The .schema
filename extension is added if needed. The schema file must reside in the library directory.
If this option is missing, the appropriate schema is either derived from the input data
set (if one is specified), or the user is asked to select a schema from those found in the
library directory.
*/
class SmileToolCommandLineHandler final
{
public:
    /** This function processes the command line arguments and invokes the appropriate high-level
        functions to perform the actions requested by the user. The function returns an appropriate
        program exit value. */
    static int perform();
};

////////////////////////////////////////////////////////////////////

#endif

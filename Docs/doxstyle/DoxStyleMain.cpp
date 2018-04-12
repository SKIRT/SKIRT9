/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DoxStyleMain.hpp"
#include "Chunk.hpp"
#include "SignalHandler.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // initialize the system
    System system(argc, argv);
    SignalHandler::InstallSignalHandlers();

    // perform a streamline operation from stdin to stdout
    Chunk styler;
    styler.readFromConsole();
    styler.streamline();
    styler.writeToConsole();
}

////////////////////////////////////////////////////////////////////

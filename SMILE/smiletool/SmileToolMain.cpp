/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SignalHandler.hpp"
#include "SmileToolCommandLineHandler.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // Initialize the system
    System system(argc, argv);
    SignalHandler::InstallSignalHandlers();

    // Handle and act on command line arguments
    return SmileToolCommandLineHandler::perform();
}

////////////////////////////////////////////////////////////////////

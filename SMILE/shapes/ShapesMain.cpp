/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SignalHandler.hpp"
#include "ShapesCommandLineHandler.hpp"
#include "ShapeRegistry.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    // Initialize the system
    System system(argc, argv);
    SignalHandler::InstallSignalHandlers();

    // Add all shapes to the item registry
    ShapeRegistry registry;

    // Handle and act on command line arguments
    return ShapesCommandLineHandler::perform();
}

////////////////////////////////////////////////////////////////////

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BuildInfo.hpp"
#include "FatalError.hpp"
#include "MainWindow.hpp"
#include "SignalHandler.hpp"
#include "System.hpp"
#include <QApplication>

////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // Construct application object
    QApplication app(argc, argv);
    app.setOrganizationName("Ghent University, Astronomical Observatory, SKIRT research team");
    app.setOrganizationDomain("skirt.ugent.be");
    app.setApplicationName("MakeUp");
    app.setApplicationVersion(QString::fromStdString(BuildInfo::projectVersion() + " (" + BuildInfo::codeVersion() + " "
                                                     + BuildInfo::timestamp() + ")"));

    // Initialize the system (do this after constructing QApplication to ensure proper locale)
    System system(argc, argv);
    SignalHandler::InstallSignalHandlers();

    // Show the main (and only) window
    MainWindow wizard;
    wizard.show();

    // Execute the event loop
    try
    {
        return app.exec();
    }
    catch (FatalError& error)
    {
        for (auto line : error.message()) System::log(line, System::LogLevel::Error);
    }
    catch (const std::exception& except)
    {
        System::log("Standard Library Exception: " + string(except.what()), System::LogLevel::Error);
    }
    return EXIT_FAILURE;
}

////////////////////////////////////////////////////////////////////

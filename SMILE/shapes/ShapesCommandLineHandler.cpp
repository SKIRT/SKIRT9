/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ShapesCommandLineHandler.hpp"
#include "BuildInfo.hpp"
#include "Canvas.hpp"
#include "CommandLineArguments.hpp"
#include "Console.hpp"
#include "FatalError.hpp"
#include "SchemaDef.hpp"
#include "ShapeCanvas.hpp"
#include "ShapeRegistry.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "XmlHierarchyCreator.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // generates an image with the specified output path for the parameter file with the specified input path;
    // returns true if successful, false otherwise
    bool doGenerateImage(string inPath, string outPath)
    {
        Console::info("Reading shapes parameter file '" + inPath + "' ...");
        Console::info("Generating shapes image file '" + outPath + "' ...");

        auto schema = ShapeRegistry::getSchemaDef();
        auto canvas = XmlHierarchyCreator::readFile(schema, inPath);
        dynamic_cast<ShapeCanvas*>(canvas.get())->paintAndSave(outPath);
        return true;
    }

    // generates a SMILE schema file with the specified file path;
    // returns true if successful, false otherwise
    bool doGenerateSchemaFile(string outPath)
    {
        Console::info("Generating SMILE schema file '" + outPath + "' ...");

        auto schema = ShapeRegistry::getSchemaDef();
        schema->save(outPath, "Shapes " + BuildInfo::projectVersion() + " " + BuildInfo::timestamp());
        return true;
    }
}

////////////////////////////////////////////////////////////////////

int ShapesCommandLineHandler::perform()
{
    // Catch and properly report any exceptions
    try
    {
        Console::warning("Welcome to 'shapes', a SMILE example " + BuildInfo::projectVersion() + " "
                         + BuildInfo::timestamp());
        Console::info("Running on host " + System::hostname() + " for " + System::username());
        Console::info("Executable path " + System::executablePath());

        // Attempt parsing both forms of command line arguments
        CommandLineArguments args1(System::arguments(), "-o*");
        CommandLineArguments args2(System::arguments(), "-x*");

        // First form: generate an image from a parameter file
        if (args1.isValid() && args1.filepaths().size() == 1)
        {
            // Get input file path; verify that the file exists; add filename extension if needed
            string inPath = args1.filepaths()[0];
            if (!System::ifstream(inPath))
            {
                inPath = StringUtils::addExtension(inPath, "shapes");
                if (!System::ifstream(inPath))
                {
                    Console::error("Can't find or open shapes parameter file '" + inPath + "'");
                    return EXIT_FAILURE;
                }
            }

            // Get or construct output file path
            string outPath = args1.value("-o");
            if (outPath.empty())
            {
                outPath = inPath;
                if (StringUtils::endsWith(StringUtils::toLower(outPath), ".shapes"))
                    outPath.erase(outPath.length() - 7);
            }
            outPath = StringUtils::addExtension(outPath, "tiff");

            if (!doGenerateImage(inPath, outPath)) return EXIT_FAILURE;
        }

        // Second form: generate a SMILE schema file
        else if (args2.isValid() && !args2.hasFilepaths() && args2.hasOptions())
        {
            string outPath = StringUtils::joinPaths(args2.value("-x"), "shapes.smile");
            if (!doGenerateSchemaFile(outPath)) return EXIT_FAILURE;
        }

        // Invalid command line arguments
        else
        {
            Console::error("Invalid command line arguments. Usage synopsis:");
            Console::warning("shapes <shapes_dataset_filepath> [-o <tiff_output_filepath>]");
            Console::warning("shapes -x <library_dirpath>");
            return EXIT_FAILURE;
        }

        // Report successful completion
        Console::success("Successful completion");

        // Report memory usage
        size_t avail = System::availableMemory();
        size_t peak = System::peakMemoryUsage();
        double peakPerCent = 100. * static_cast<double>(peak) / static_cast<double>(avail);
        Console::info("Peak memory usage: " + StringUtils::toMemSizeString(peak) + " ("
                      + StringUtils::toString(peakPerCent, 'f', 1) + "% of " + StringUtils::toMemSizeString(avail)
                      + ")");
        return EXIT_SUCCESS;
    }
    catch (const FatalError& error)
    {
        for (auto line : error.message()) Console::error(line);
    }
    catch (const std::exception& except)
    {
        Console::error("Standard Library Exception: " + string(except.what()));
    }
    return EXIT_FAILURE;
}

////////////////////////////////////////////////////////////////////

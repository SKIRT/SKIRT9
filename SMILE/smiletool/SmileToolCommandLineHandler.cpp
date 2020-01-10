/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SmileToolCommandLineHandler.hpp"
#include "BuildInfo.hpp"
#include "CommandLineArguments.hpp"
#include "Console.hpp"
#include "ConsoleHierarchyCreator.hpp"
#include "FatalError.hpp"
#include "Item.hpp"
#include "LatexHierarchyWriter.hpp"
#include "SchemaDef.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "XmlHierarchyCreator.hpp"
#include "XmlHierarchyWriter.hpp"

////////////////////////////////////////////////////////////////////

int SmileToolCommandLineHandler::perform()
{
    // Catch and properly report any exceptions
    try
    {
        Console::warning("Welcome to the SMILE console tool " + BuildInfo::projectVersion() + " "
                         + BuildInfo::timestamp());

        // Process and validate the command line arguments
        CommandLineArguments args(System::arguments(), "-i* -o* -t* -l* -s*");
        if (!args.isValid() || args.hasFilepaths())
        {
            Console::error("Invalid command line arguments. Usage synopsis:");
            Console::warning("smiletool [-i <input_dataset_filepath>]");
            Console::warning("          [-o <output_dataset_filepath>] [-t <output_latex_filepath>]");
            Console::warning("          [-l <library_dirpath>] [-s <schema_filename>]");
            return EXIT_FAILURE;
        }

        // Make a list of schema files residing in the library directory, and a list of corresponding titles
        string libraryPath = args.value("-l");
        vector<string> schemaNames;
        vector<string> schemaTitles;
        for (string name : System::filesInDirectory(libraryPath))
        {
            if (StringUtils::endsWith(name, ".smile"))
            {
                string title = SchemaDef::getSchemaTitle(StringUtils::joinPaths(libraryPath, name));
                if (!title.empty())
                {
                    schemaNames.push_back(name);
                    schemaTitles.push_back(title);
                }
            }
        }
        if (schemaTitles.empty())
        {
            Console::error("There are no SMILE schema files in the library directory '" + libraryPath + "'");
            return EXIT_FAILURE;
        }

        // Obtain input dataset filepath from the command line, if specified
        string infilename = args.value("-i");

        // Obtain SMILE schema filename from the command line, from the input dataset, or by prompting the user
        string schemaName = StringUtils::addExtension(args.value("-s"), "smile");
        if (!schemaName.empty())
        {
            if (StringUtils::indexOf(schemaNames, schemaName) < 0)
            {
                Console::error("The specified SMILE schema file '" + schemaName
                               + "' is not found in the library directory '" + libraryPath + "'");
                return EXIT_FAILURE;
            }
        }
        else if (!infilename.empty())  // loading existing dataset
        {
            // search for compatible schema in the libary
            for (const auto& name : schemaNames)
            {
                if (SchemaDef::isCompatible(StringUtils::joinPaths(libraryPath, name), infilename))
                {
                    schemaName = name;
                    break;
                }
            }

            // if not found, produce an error
            if (schemaName.empty())
            {
                Console::error("The SMILE schema library does not contain a schema compatible with input dataset '"
                               + infilename + "'");
                return EXIT_FAILURE;
            }
        }
        else  // constructing new dataset
        {
            auto schemaIndex = Console::promptForChoice("the file type to be created", schemaTitles, false, 0, true,
                                                        "or zero to exit");
            if (schemaIndex < 0)
            {
                Console::warning("Session aborted by user");
                return EXIT_FAILURE;
            }
            schemaName = schemaNames[schemaIndex];
        }

        // Load the selected SMILE schema
        SchemaDef schema(StringUtils::joinPaths(libraryPath, schemaName));

        // Verify that the input dataset exists (if it is specified).
        if (!infilename.empty() && !System::ifstream(infilename))
        {
            infilename = StringUtils::addExtension(infilename, schema.schemaExtension());
            if (!System::ifstream(infilename))
            {
                Console::error("Can't find or open dataset '" + infilename + "'");
                return EXIT_FAILURE;
            }
        }

        // Inform the user of the proceedings
        if (!infilename.empty())
            Console::warning("Loading " + infilename + " (using " + schemaName + ")");
        else
            Console::warning("Creating " + schema.schemaTitle() + " (using " + schemaName + ")");
        if (!schema.schemaProducer().empty())
            Console::info(schema.schemaName() + " schema produced by " + schema.schemaProducer());

        // Load a SMILE dataset from file or create a new one through a console Q&A session
        auto top = !infilename.empty() ? XmlHierarchyCreator::readFile(&schema, infilename)
                                       : ConsoleHierarchyCreator::create(&schema);

        // Obtain an output file name for the dataset file from the command line, or by prompting the user
        string outfilename = StringUtils::addExtension(args.value("-o"), schema.schemaExtension());
        if (outfilename.empty() && (!args.isPresent("-i") || !args.isPresent("-t")))
        {
            while (true)
            {
                // get a file name
                outfilename = Console::promptForString(
                    "Enter the name of the " + schema.schemaExtension() + " file to be saved", false, string());
                outfilename = StringUtils::addExtension(outfilename, schema.schemaExtension());
                // verify that the file does not yet exists
                // (we test whether the file can be opened, which is the best we can do in standard C++14)
                if (System::ifstream(outfilename))
                    Console::error("File '" + outfilename + "' already exists; enter another name");
                else
                    break;
            }
        }

        // If requested, save the dataset as XML and inform the user
        if (!outfilename.empty())
        {
            XmlHierarchyWriter::write(top.get(), &schema, outfilename,
                                      "smiletool " + BuildInfo::projectVersion() + " " + BuildInfo::timestamp());
            Console::success("Saved " + schema.schemaTitle() + " to '" + outfilename + "'");
        }

        // If requested, write the dataset as LaTeX and inform the user
        if (args.isPresent("-t"))
        {
            string latexfilename = StringUtils::addExtension(args.value("-t"), "tex");
            string datasetname = infilename;
            if (datasetname.empty()) datasetname = outfilename;
            if (datasetname.empty()) datasetname = "SMILE";

            LatexHierarchyWriter::write(top.get(), &schema, latexfilename, datasetname);
            Console::success("Saved LaTeX overview to '" + latexfilename + "'");
        }

        // Report successful completion
        Console::success("Successful completion");
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

/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustGrainPopulationsProbe.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "MediumSystem.hpp"
#include "MultiGrainPopulationInterface.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // return the number of digits in a number in the range [0, 99999]
    int numDigits(int number)
    {
        if (number < 10) return 1;
        if (number < 100) return 2;
        if (number < 1000) return 3;
        if (number < 10000) return 4;
        return 5;
    }
}

////////////////////////////////////////////////////////////////////

void DustGrainPopulationsProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        auto units = find<Units>();

        // locate the medium system
        auto ms = find<MediumSystem>();
        int numMedia = ms->numMedia();

        // create a seperate file for each medium
        for (int h = 0; h != numMedia; ++h)
        {
            // get the mix and skip mixes that don't offer multiple dust grain populations
            auto mix = dynamic_cast<const MultiGrainPopulationInterface*>(ms->media()[h]->mix());
            if (mix)
            {
                // create a text file
                TextOutFile out(this, itemName() + "_grainpops_" + std::to_string(h), "grain populations");

                // write the header
                double totalmass = mix->totalMass();
                out.writeLine("# Medium component " + std::to_string(h) + " -- dust mass per hydrogen atom: "
                              + StringUtils::toString(units->obulkmass(totalmass), 'e', 9) + " " + units->ubulkmass());
                out.addColumn("grain population index");
                out.addColumn("grain material type", "-");
                out.addColumn("dust mass as a percentage of total", "%");
                out.addColumn("dust mass per hydrogen atom", units->ubulkmass());
                out.addColumn("dust mass per hydrogen mass");
                out.addColumn("minimum grain size", units->ugrainsize());
                out.addColumn("maximum grain size", units->ugrainsize());

                // determine the number of digits in the largest population index
                int numPops = mix->numPopulations();
                int indexWidth = numDigits(numPops - 1);

                // determine the number of characters in the longest material type description
                size_t typeWidth = 1;
                for (int c = 0; c != numPops; ++c) typeWidth = max(typeWidth, mix->populationGrainType(c).size());

                // write the columns
                for (int c = 0; c != numPops; ++c)
                {
                    string line = StringUtils::toString(c, 'd', 0, indexWidth);
                    line += " " + StringUtils::padRight(mix->populationGrainType(c), typeWidth);

                    double mass = mix->populationMass(c);
                    line += " " + StringUtils::toString(100. * mass / totalmass, 'f', 4, 8);
                    line += " " + StringUtils::toString(units->obulkmass(mass), 'e', 9);
                    line += " " + StringUtils::toString(mass / Constants::Mproton(), 'e', 9);

                    Range range = mix->populationSizeRange(c);
                    line += " " + StringUtils::toString(units->ograinsize(range.min()), 'e', 9);
                    line += " " + StringUtils::toString(units->ograinsize(range.max()), 'e', 9);

                    out.writeLine(line);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

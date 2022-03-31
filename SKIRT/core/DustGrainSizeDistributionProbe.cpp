/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustGrainSizeDistributionProbe.hpp"
#include "Configuration.hpp"
#include "GrainSizeDistribution.hpp"
#include "MediumSystem.hpp"
#include "MultiGrainPopulationInterface.hpp"
#include "NR.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void DustGrainSizeDistributionProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        auto units = find<Units>();

        // locate the medium system
        auto ms = find<MediumSystem>();
        int numMedia = ms->numMedia();

        // loop over each medium and each population
        // skipping mixes that don't offer multiple dust grain populations
        for (int h = 0; h != numMedia; ++h)
        {
            auto mix = dynamic_cast<const MultiGrainPopulationInterface*>(ms->media()[h]->mix());
            if (mix)
            {
                int numPops = mix->numPopulations();
                for (int c = 0; c != numPops; ++c)
                {
                    // get the size distribution
                    auto sd = mix->populationSizeDistribution(c);

                    // create a text file
                    TextOutFile out(this, itemName() + "_grainsizes_" + std::to_string(h) + "_" + std::to_string(c),
                                    "grain size distribution");

                    // write the header
                    out.writeLine("# Dust grain size distribution");
                    out.addColumn("grain size", units->ugrainsize());
                    out.addColumn("size distribution", units->upergrainsize());

                    // construct a grain size grid
                    Array av;
                    NR::buildLogGrid(av, sd->amin(), sd->amax(), _numSamples - 1);

                    // write the columns
                    for (int i = 0; i != _numSamples; ++i)
                    {
                        out.writeRow(units->ograinsize(av[i]), units->opergrainsize(sd->dnda(av[i])));
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

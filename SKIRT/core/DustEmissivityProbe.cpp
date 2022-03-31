/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustEmissivityProbe.hpp"
#include "ArrayTable.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Indices.hpp"
#include "InstrumentWavelengthGridProbe.hpp"
#include "MediumSystem.hpp"
#include "PlanckFunction.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // This function returns the Mathis field J_lambda (Mathis et al. 1983)
    // discretized on the simulation's radiation field wavelength grid.
    // Note that the Mathis recipe describes the field as 4pi J_lambda,
    // whereas this function returns the mean intensity J_lambda per steradian.
    Array mathis(Probe* probe)
    {
        auto wavelengthGrid = probe->find<Configuration>()->radiationFieldWLG();
        int numWavelengths = wavelengthGrid->numBins();
        int ellA = max(0, wavelengthGrid->bin(0.0912e-6));
        int ellB = max(ellA, wavelengthGrid->bin(0.110e-6));
        int ellC = max(ellB, wavelengthGrid->bin(0.134e-6));
        int ellD = max(ellC, wavelengthGrid->bin(0.250e-6));

        const double Wv[] = {1e-14, 1e-13, 4e-13};
        const double Tv[] = {7500, 4000, 3000};

        Array Jv(numWavelengths);
        for (int ell = ellA + 1; ell <= ellB; ell++)
            Jv[ell] = 3069. * pow(wavelengthGrid->wavelength(ell) * 1e6, 3.4172);
        for (int ell = ellB + 1; ell <= ellC; ell++) Jv[ell] = 1.627;
        for (int ell = ellC + 1; ell <= ellD; ell++)
            Jv[ell] = 0.0566 * pow(wavelengthGrid->wavelength(ell) * 1e6, -1.6678);
        for (int i = 0; i < 3; i++)
        {
            PlanckFunction B(Tv[i]);
            for (int ell = ellD + 1; ell < numWavelengths; ell++) Jv[ell] += Wv[i] * B(wavelengthGrid->wavelength(ell));
        }
        return Jv;
    }

    // This function returns the blackbody field B_lambda(T)
    // discretized on the simulation's radiation field wavelength grid.
    Array blackbody(Probe* probe, double T)
    {
        auto wavelengthGrid = probe->find<Configuration>()->radiationFieldWLG();
        int numWavelengths = wavelengthGrid->numBins();

        PlanckFunction B(T);
        Array Jv(numWavelengths);
        for (int ell = 0; ell < numWavelengths; ell++) Jv[ell] = B(wavelengthGrid->wavelength(ell));
        return Jv;
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // this function writes one of the output files for this probe;
    // Jv must be discretized on the simulation's radiation field wavelength grid
    void writeEmissivitiesForField(Probe* probe, const Array& Jv, string name, string title)
    {
        auto ms = probe->find<MediumSystem>();
        auto units = probe->find<Units>();
        auto wavelengthGrid = probe->find<Configuration>()->dustEmissionWLG();

        // construct a list of indices and material mixes for medium components that actually contain dust
        vector<int> hv;
        vector<const MaterialMix*> mixv;
        for (int h = 0; h != ms->numMedia(); ++h)
            if (ms->isDust(h))
            {
                hv.push_back(h);
                mixv.push_back(ms->media()[h]->mix());
            }

        // calculate the emissivity for each representative dust mix
        ArrayTable<2> evv(hv.size(), 0);
        for (size_t i = 0; i != hv.size(); ++i) evv(i) = mixv[i]->emissivity(Jv);

        // create an output text file
        TextOutFile file(probe, probe->itemName() + "_" + name, "dust emissivities for " + title);

        // write the header
        file.writeLine("# Dust emissivities for input field " + title);
        file.addColumn("wavelength; " + units->swavelength(), units->uwavelength());
        for (int h : hv) file.addColumn("lambda*j_lambda for dust in medium component " + std::to_string(h), "W/sr/H");

        // write the emissivity for each dust mix to file
        for (int ell : Indices(wavelengthGrid->numBins(), units->rwavelength()))
        {
            double lambda = wavelengthGrid->wavelength(ell);
            vector<double> values({units->owavelength(lambda)});
            for (size_t i = 0; i != hv.size(); ++i) values.push_back(lambda * evv(i, ell + 1));
            // (add 1 to ell to skip leftmost wavelength grid border point included in evv)
            file.writeRow(values);
        }
    }
}

////////////////////////////////////////////////////////////////////

void DustEmissivityProbe::probe()
{
    if (find<Configuration>()->hasDustEmission() && find<MediumSystem>()->hasDust())
    {
        // write emissivities for a range of scaled Mathis ISRF input fields
        {
            Array Jv = mathis(this);
            for (int i = -4; i < 7; i++)
            {
                double U = pow(10., i);
                writeEmissivitiesForField(this, U * Jv, "Mathis_U_" + StringUtils::toString(U, 'e', 0),
                                          StringUtils::toString(U, 'g') + " * Mathis ISRF");
            }
        }

        // write emissivities for a range of diluted black body input fields
        {
            const int Tv[] = {3000, 6000, 9000, 12000, 15000, 18000};
            const double Dv[] = {8.28e-12, 2.23e-13, 2.99e-14, 7.23e-15, 2.36e-15, 9.42e-16};
            for (int i = 0; i < 6; i++)
            {
                writeEmissivitiesForField(
                    this, Dv[i] * blackbody(this, Tv[i]), "BlackBody_T_" + StringUtils::toString(Tv[i], 'd', 0, 5, '0'),
                    StringUtils::toString(Dv[i], 'e', 2) + " * B(" + StringUtils::toString(Tv[i]) + "K)");
            }
        }

        // if requested, also output the wavelength grid
        if (writeWavelengthGrid())
        {
            InstrumentWavelengthGridProbe::writeWavelengthGrid(this, find<Configuration>()->dustEmissionWLG(),
                                                               itemName() + "_wavelengths",
                                                               "emission spectrum wavelengths");
        }
    }
}

////////////////////////////////////////////////////////////////////

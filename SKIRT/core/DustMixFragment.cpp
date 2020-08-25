/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustMixFragment.hpp"

//////////////////////////////////////////////////////////////////////

DustMixFragment::DustMixFragment(SimulationItem* parent, ScatteringMode scatteringMode,
                                 const GrainPopulation* population)
    : _scatteringMode{scatteringMode}, _population{population}
{
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

namespace
{
    // This class decorates an arbitrary size distribution so that it becomes limited to the given size range.
    class SizeDistributionFragment : public GrainSizeDistribution
    {
    public:
        explicit SizeDistributionFragment(SimulationItem* parent, const GrainSizeDistribution* original, double amin,
                                          double amax)
            : _original{original}, _amin{amin}, _amax{amax}
        {
            parent->addChild(this);
            setup();
        }
        double amin() const override { return _amin; }
        double amax() const override { return _amax; }
        double dnda(double a) const override { return _original->dnda(a); }

    private:
        const GrainSizeDistribution* _original;
        double _amin;
        double _amax;
    };
}

//////////////////////////////////////////////////////////////////////

DustMixFragment::DustMixFragment(SimulationItem* parent, ScatteringMode scatteringMode,
                                 const GrainPopulation* population, double amin, double amax)
    : _scatteringMode{scatteringMode}
{
    // TO DO: fix normalization for new population
    _population = new GrainPopulation(this, population->composition(),
                                      new SizeDistributionFragment(this, population->sizeDistribution(), amin, amax), 1,
                                      GrainPopulation::NormalizationType::FactorOnSizeDistribution, 1.);
    parent->addChild(this);
    setup();
}

//////////////////////////////////////////////////////////////////////

void DustMixFragment::setupSelfBefore()
{
    MultiGrainDustMix::setupSelfBefore();

    addPopulation(_population);
}

//////////////////////////////////////////////////////////////////////

DustMix::ScatteringMode DustMixFragment::scatteringMode() const
{
    return _scatteringMode;
}

////////////////////////////////////////////////////////////////////

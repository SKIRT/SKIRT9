/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustMixFragment.hpp"
#include "StringUtils.hpp"

//////////////////////////////////////////////////////////////////////

string DustMixFragment::type() const
{
    return "FragmentDustMixDecorator";
}

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
                                 const GrainPopulation* population, double amin, double amax, double normalization)
    : _scatteringMode{scatteringMode}
{
    _population = new GrainPopulation(this, population->composition(),
                                      new SizeDistributionFragment(this, population->sizeDistribution(), amin, amax), 1,
                                      GrainPopulation::NormalizationType::FactorOnSizeDistribution, normalization);
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

void DustMixFragment::setupSelfAfter()
{
    MultiGrainDustMix::setupSelfAfter();

    // graphite or silicate?
    string name = _population->composition()->name();
    _isGraphite =
        StringUtils::contains(name, "Gra") || StringUtils::contains(name, "PAH") || StringUtils::contains(name, "CM20");

    // integrate over the grain size distribution to calculate the average grain radius
    double logamin = log10(_population->sizeDistribution()->amin());
    double logamax = log10(_population->sizeDistribution()->amax());
    int numSizes = max(3., 100. * (logamax - logamin));
    double dloga = (logamax - logamin) / (numSizes - 1);
    double sum1 = 0.;
    double sum2 = 0.;
    for (int i = 0; i != numSizes; ++i)
    {
        double w = (i == 0 || i == numSizes - 1) ? 0.5 : 1.;
        double a = pow(10, logamin + i * dloga);
        double da = a * M_LN10 * dloga;
        double dnda = _population->sizeDistribution()->dnda(a);
        sum1 += w * dnda * a * da;
        sum2 += w * dnda * da;
    }
    _grainRadius = sum2 ? sum1 / sum2 : 0.;
}

//////////////////////////////////////////////////////////////////////

DustMix::ScatteringMode DustMixFragment::scatteringMode() const
{
    return _scatteringMode;
}

//////////////////////////////////////////////////////////////////////

bool DustMixFragment::isGraphite() const
{
    return _isGraphite;
}

//////////////////////////////////////////////////////////////////////

double DustMixFragment::grainRadius() const
{
    return _grainRadius;
}

////////////////////////////////////////////////////////////////////

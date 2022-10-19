/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEMPERATUREPROBE_HPP
#define TEMPERATUREPROBE_HPP

#include "SpatialGridWhenFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** TemperatureProbe probes the indicative temperature of the medium as discretized on the spatial
    grid of the simulation. The meaning of the indicative temperature depends heavily on the
    material type; see the discussion below. The probe can be used with any Form subclass. When
    associated with a form that projects the quantity along a path, the indicative temperature
    value is density-weighted.

    The user can select the aggregation level, i.e. whether to produce an output file per medium
    component or per medium type (dust, electrons, gas). There is also an option to decide whether
    the probe should be performed after setup or after the full simulation run. See the discussion
    below for more information on these options.

    <b>Indicative temperature for a single medium component</b>

    The definition of the indicative temperature \f${\bar{T}}\f$ depends heavily on the material
    type. Note that the indicative temperature almost never corresponds to a physical temperature;
    it is intended to provide an indication, as implied by its name.

    For a dust component, the indicative temperature in a given spatial cell is obtained by solving
    the energy balance equation under local thermal equilibrium (LTE) assumptions for a single
    representative grain of the associated dust mix (even if the dust mix includes multiple grain
    composition types and/or grains size bins). This implies that the indicative temperature for a
    dust component can be calculated only after the radiation field has been calculated by the
    simulation.

    For an electron or gas component, the indicative temperature is equal to the imported or
    default temperature provided as part of the input model and stored in the medium state for each
    spatial cell during setup. In principle, a gas mix may update the temperature during the
    simulation, however none of the current or foreseeable gas mixes do this.

    In all cases, if a cell does not contain any material for the requested component, the
    indicative temperature is taken to be zero.

    <b>Indicative temperature for all medium components of a given material type</b>

    The aggregated indicative temperature for all medium components of a given material type (dust,
    electrons, gas) is obtained by averaging the individual indicative temperatures, as defined
    above, weighted by the relative masses in the cell under consideration. In formula form, for a
    spatial cell \f$m\f$ with components \f$h\f$ of the specified material type, the indicative
    dust temperature is defined as \f[{\bar{T}}_m = \frac{\sum_h \rho_{m,h}\,{\bar{T}}_{m,h}}
    {\sum_h \rho_{m,h}} \f] where \f${\bar{T}}_{m,h}\f$ is the indicative temperature for each
    component \f$h\f$ and \f$\rho_{m,h}\f$ is the corresponding (mass or number) density.

    <b>Indicative temperature along a path</b>

    When a TemperatureProbe is associated with a Form subclass that projects the probed quantity
    along a path, the indicative temperature is similarly density-weighted over the path, i.e.
    \f[{\bar{T}} = \frac{\sum_m \rho_m \,{\bar{T}}_m \,\Delta s_m} {\sum_m \rho_m \,\Delta s_m}\f]
    where \f$\Delta s_m\f$ is the length of the path segment in each crossed cell \f$m\f$.

    <b>Probing after setup or at the end of the simulation run</b>

    As mentioned above, the indicative dust temperature can be obtained only after the simulation
    has calculated the radiation field. The probe will skip any dust components if it is invoked
    after setup or if the simulation does not record a panchromatic radiation field.

    This consideration does not hold for electron and gas components; the indicative electron or
    gas temperature can be probed either after setup or at the end of the simulation run. */
class TemperatureProbe : public SpatialGridWhenFormProbe
{
    /** The enumeration type indicating how to aggregate the output: per medium component or per
        medium type (dust, electrons, gas). */
    ENUM_DEF(Aggregation, Component, Type)
        ENUM_VAL(Aggregation, Component, "per medium component")
        ENUM_VAL(Aggregation, Type, "per medium type (dust, electrons, gas)")
    ENUM_END()

    ITEM_CONCRETE(TemperatureProbe, SpatialGridWhenFormProbe,
                  "internal spatial grid: indicative temperature of the medium")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TemperatureProbe,
                                    "Medium&SpatialGrid&((DustMix&RadiationField&Panchromatic)|ElectronMix|GasMix)")

        PROPERTY_ENUM(aggregation, Aggregation, "how to aggregate the indicative temperature")
        ATTRIBUTE_DEFAULT_VALUE(aggregation, "Type")
        ATTRIBUTE_DISPLAYED_IF(aggregation, "Level2")
        ATTRIBUTE_INSERT(aggregation, "thisIsTemperatureProbe")  // used to set default in SpatialGridWhenFormProbe

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif

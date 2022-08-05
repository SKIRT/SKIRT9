/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSYSTEM_HPP
#define MEDIUMSYSTEM_HPP

#include "Array.hpp"
#include "DustEmissionOptions.hpp"
#include "DynamicStateOptions.hpp"
#include "IterationOptions.hpp"
#include "LyaOptions.hpp"
#include "MaterialMix.hpp"
#include "Medium.hpp"
#include "MediumState.hpp"
#include "PhotonPacketOptions.hpp"
#include "RadiationFieldOptions.hpp"
#include "SamplingOptions.hpp"
#include "SecondaryEmissionOptions.hpp"
#include "SimulationItem.hpp"
#include "SpatialGrid.hpp"
#include "Table.hpp"
class Configuration;
class PhotonPacket;
class Random;
class ShortArray;
class SpatialGridPath;
class WavelengthGrid;

//////////////////////////////////////////////////////////////////////

/** An instance of the MediumSystem class represents a complete medium system, which is the
    superposition of one or more transfer media. Each individual medium represents a spatial
    density distribution and defines the material properties of the medium at each location. While
    the specific material properties may vary with location, the fundamental material type must be
    the same throughout the spatial domain for each medium.

    In addition to the media input model, the MediumSystem class includes the spatial grid that
    tessellates the spatial domain of the simulation into cells, and manages the medium state and
    the radiation field for each spatial cell in this grid. The class therefore plays a central
    role in the simulation and it offers quite a few different type of functions.

    <b>Overall medium configuration</b>

    These functions offer basic information on the spatial grid, the medium components, and the
    material types in the medium system.

    <b>%Medium state</b>

    The medium state for each spatial cell in the simulation includes a set of \em common state
    variables shared by all medium components and a set of \em specific state variables for each
    individual medium component. See the MediumState class for more information. The functions in
    this section allow access to the common state variables for a given spatial cell.

    <b>Low-level optical properties</b>

    These functions allow retrieving absorption, scattering and extinction cross sections for a
    given spatial cell and material type as a function of wavelength, assuming fixed, predefined
    values for any quantities other than wavelength (e.g., a default temperature, no polarization,
    no kinematics). The values returned by these low-level functions may be used only during setup
    and for probing.

    <b>High-level photon life cycle</b>

    These functions support the photon life cycle by offering a generic interface for operations
    including tracing photon paths through the medium and performing scattering events. The
    functions calculate opacities and other medium properties in a given spatial cell by passing
    the incoming photon packet and the full medium state to the appropriate material mix for each
    medium component. This allows proper treatment of polarization and kinematics and supports
    dependencies on temperature or custom special state variables.

    <b>Radiation field</b>

    These functions allow storing the radiation field during the photon life cycle and retrieving
    the results after a set of photon's have been processed. The contribution to the radation field
    for each spatial cell and for each wavelength in the simulation's radiation field wavelength
    grid is traced separately for primary and secondary sources. This avoids the need for repeating
    primary emission during dust-temperature convergence iterations. At all times, the sum of the
    primary and secondary contributions represents the radiation field to be used as input for
    calculations. There is a third, temporary table that serves as a target for storing the
    secondary radiation field so that the "stable" primary and secondary tables remain available
    for calculating secondary emission spectra while shooting secondary photons through the grid.

    <b>Indicative temperature</b>

    These functions determine an indicative temperature in a given spatial cell and for a given
    material type, averaged over medium components if applicable. Depending on the material type,
    the indicative temperature may be based on the radiation field or it may be derived from a
    value given in the input model. In any case, it does not reflect a physical quantity and it
    should be used only for setup and probing purposes.

    <b>Emission</b>

    These functions support the secondary source system by offering a generic interface for
    calculating the secondary emission properties in a given spatial cell and for a given material
    type, including luminosities, spectra and polarization. */
class MediumSystem : public SimulationItem
{
    ITEM_CONCRETE(MediumSystem, SimulationItem, "a medium system")
        ATTRIBUTE_TYPE_ALLOWED_IF(MediumSystem, "!NoMedium")

        PROPERTY_ITEM(photonPacketOptions, PhotonPacketOptions, "the photon packet options")
        ATTRIBUTE_DEFAULT_VALUE(photonPacketOptions, "PhotonPacketOptions")
        ATTRIBUTE_RELEVANT_IF(media, "!NoMedium")

        PROPERTY_ITEM(lyaOptions, LyaOptions, "the Lyman-alpha line transfer options")
        ATTRIBUTE_DEFAULT_VALUE(lyaOptions, "LyaOptions")
        ATTRIBUTE_RELEVANT_IF(lyaOptions, "Lya")

        PROPERTY_ITEM(dynamicStateOptions, DynamicStateOptions, "the dynamic medium state options")
        ATTRIBUTE_DEFAULT_VALUE(dynamicStateOptions, "DynamicStateOptions")
        ATTRIBUTE_RELEVANT_IF(dynamicStateOptions, "IteratePrimary")

        PROPERTY_ITEM(radiationFieldOptions, RadiationFieldOptions, "the radiation field options")
        ATTRIBUTE_DEFAULT_VALUE(radiationFieldOptions, "RadiationFieldOptions")
        ATTRIBUTE_RELEVANT_IF(radiationFieldOptions, "!NoMedium")

        PROPERTY_ITEM(secondaryEmissionOptions, SecondaryEmissionOptions, "the secondary emission options")
        ATTRIBUTE_DEFAULT_VALUE(secondaryEmissionOptions, "SecondaryEmissionOptions")
        ATTRIBUTE_RELEVANT_IF(secondaryEmissionOptions, "Emission")

        PROPERTY_ITEM(iterationOptions, IterationOptions, "the primary and/or secondary emission iteration options")
        ATTRIBUTE_DEFAULT_VALUE(iterationOptions, "IterationOptions")
        ATTRIBUTE_RELEVANT_IF(iterationOptions, "IteratePrimary|IterateSecondary")

        PROPERTY_ITEM(dustEmissionOptions, DustEmissionOptions, "the dust emission options")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissionOptions, "DustEmissionOptions")
        ATTRIBUTE_RELEVANT_IF(dustEmissionOptions, "DustEmission")

        PROPERTY_ITEM_LIST(media, Medium, "the transfer media")
        ATTRIBUTE_DEFAULT_VALUE(media, "GeometricMedium")
        ATTRIBUTE_REQUIRED_IF(media, "!NoMedium")

        PROPERTY_ITEM(samplingOptions, SamplingOptions, "the spatial grid sampling options")
        ATTRIBUTE_DEFAULT_VALUE(samplingOptions, "SamplingOptions")

        PROPERTY_ITEM(grid, SpatialGrid, "the spatial grid")
        ATTRIBUTE_DEFAULT_VALUE(grid,
                                "Dimension3:PolicyTreeSpatialGrid;Dimension2:Cylinder2DSpatialGrid;Sphere1DSpatialGrid")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates and stores initial state information for each spatial cell,
        including the cell volume and the number density for each medium as defined by the input
        model. If needed for the simulation's configuration, it also allocates one or two radiation
        field data tables that have a bin for each spatial cell in the simulation and for each bin
        in the wavelength grid returned by the Configuration::radiationFieldWLG() function. */
    void setupSelfAfter() override;

    //=============== Overall medium configuration ===================

public:
    /** This function returns the dimension of the medium system, which depends on the (lack of)
        symmetry in the geometries of the media it contains (\em not including the spatial grid). A
        value of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of these
        symmetries. The medium with the least symmetry (i.e. the highest dimension) determines the
        result for the whole system. */
    int dimension() const;

    /** This function returns the dimension of the spatial grid held by the medium system. A value
        of 1 means spherical symmetry, 2 means axial symmetry and 3 means none of these symmetries.
        */
    int gridDimension() const;

    /** This function returns the number of media in the medium system. The returned value is valid
        only after setup has been performed. */
    int numMedia() const;

    /** This function returns the number of cells in the spatial grid held by the medium system.
        The returned value is valid only after setup has been performed. */
    int numCells() const;

    /** This function returns the material mix corresponding to the medium component with index
        \f$h\f$ in spatial cell with index \f$m\f$. */
    const MaterialMix* mix(int m, int h) const;

    /** This function returns true if at least one of the media in the medium system has the
        specified fundamental material type (i.e. dust, electrons, or gas). */
    bool hasMaterialType(MaterialMix::MaterialType type) const;

    /** This function returns true if at least one of the media in the medium system contains dust.
        */
    bool hasDust() const { return hasMaterialType(MaterialMix::MaterialType::Dust); }

    /** This function returns true if at least one of the media in the medium system contains
        electrons. */
    bool hasElectrons() const { return hasMaterialType(MaterialMix::MaterialType::Electrons); }

    /** This function returns true if at least one of the media in the medium system contains gas.
        */
    bool hasGas() const { return hasMaterialType(MaterialMix::MaterialType::Gas); }

    /** This function returns true if the medium component with index \f$h\f$ has the specified
        fundamental material type (i.e. dust, electrons, or gas). */
    bool isMaterialType(MaterialMix::MaterialType type, int h) const;

    /** This function returns true if the medium component with index \f$h\f$ contains dust. */
    bool isDust(int h) const { return isMaterialType(MaterialMix::MaterialType::Dust, h); }

    /** This function returns true if the medium component with index \f$h\f$ contains electrons.
        */
    bool isElectrons(int h) const { return isMaterialType(MaterialMix::MaterialType::Electrons, h); }

    /** This function returns true if the medium component with index \f$h\f$ contains gas. */
    bool isGas(int h) const { return isMaterialType(MaterialMix::MaterialType::Gas, h); }

    /** This function returns a list of indices \f$h\f$ for media components that contain dust. */
    const vector<int>& dustMediumIndices() const { return _dust_hv; }

    /** This function returns a list of indices \f$h\f$ for media components that contain gas. */
    const vector<int>& gasMediumIndices() const { return _gas_hv; }

    /** This function returns a list of indices \f$h\f$ for media components that contain electrons. */
    const vector<int>& electronMediumIndices() const { return _elec_hv; }

    //=============== Input model ===================

public:
    /** This function returns the total mass density of all dust medium components at the specified
        position in the input model (i.e. directly querying the medium components rather than the
        gridded medium state). */
    double dustMassDensity(Position bfr) const;

    /** This function returns the total number density of all electron medium components at the specified
        position in the input model (i.e. directly querying the medium components rather than the
        gridded medium state). */
    double electronNumberDensity(Position bfr) const;

    /** This function returns the total number density of all gas medium components at the specified
        position in the input model (i.e. directly querying the medium components rather than the
        gridded medium state). */
    double gasNumberDensity(Position bfr) const;

    //=============== Medium state ===================

public:
    /** This function returns the volume of the spatial cell with index \f$m\f$. */
    double volume(int m) const;

    /** This function returns the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in
        spatial cell with index \f$m\f$. If there are multiple media components, the aggregate bulk
        velocity \f${\boldsymbol{v}}\f$ is determined by averaging the respective bulk velocities
        over the corresponding number densities, \f[{\boldsymbol{v}} = \frac{\sum_h n_h
        {\boldsymbol{v}}_h} {\sum_h n_h}.\f] If no medium component specifies a bulk velocity, this
        function returns the null vector. */
    Vec bulkVelocity(int m) const;

    /** This function returns the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell with
        index \f$m\f$. At most one medium component is allowed to specify a magnetic field. If no
        medium component specifies a magnetic field, this function returns the null vector. */
    Vec magneticField(int m) const;

    /** This function returns the total mass density of all medium components in spatial cell with
        index \f$m\f$. */
    double massDensity(int m) const;

    /** This function returns the total mass density of all dust medium components in spatial cell
        with index \f$m\f$. */
    double dustMassDensity(int m) const;

    /** This function returns the total number density of all electron medium components in spatial
        cell with index \f$m\f$. */
    double electronNumberDensity(int m) const;

    /** This function returns the total number density of all gas medium components in spatial cell
        with index \f$m\f$. */
    double gasNumberDensity(int m) const;

    /** This function returns the number density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double numberDensity(int m, int h) const;

    /** This function returns the mass density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double massDensity(int m, int h) const;

    /** This function returns the metallicity \f$Z\f$ of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. If the specified medium component does not have the
        metallicity specific state variable, the behavior of this function is undefined. */
    double metallicity(int m, int h) const;

    /** This function returns the temperature \f$T\f$ of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. If the specified medium component does not have the
        temperature specific state variable, the behavior of this function is undefined. */
    double temperature(int m, int h) const;

    /** This function returns the value of the custom specific state variable with index \f$i\f$ of
        the medium component with index \f$h\f$ in the spatial cell with index \f$m\f$. If the
        specified medium component does not have a custom variable with the specified index, the
        behavior of this function is undefined. */
    double custom(int m, int h, int i) const;

    //=============== Low-level optical properties ===================

public:
    /** This function returns the absorption opacity \f$k_h^\text{abs}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacityAbs(double lambda, int m, int h) const;

    /** This function returns the scattering opacity \f$k_h^\text{sca}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacitySca(double lambda, int m, int h) const;

    /** This function returns the extinction opacity \f$k_h^\text{ext}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacityExt(double lambda, int m, int h) const;

    /** This function returns the absorption opacity \f$k^\text{abs}=\sum_h k_h^\text{abs}\f$
        summed over all medium components with the specified material type at wavelength
        \f$\lambda\f$ in spatial cell with index \f$m\f$. Because no photon packet is provided,
        default values are used for any relevant incoming photon packet properties. For example,
        the radiation is assumed to be unpolarized. */
    double opacityAbs(double lambda, int m, MaterialMix::MaterialType type) const;

    /** This function returns the extinction opacity \f$k^\text{ext}=\sum_h k_h^\text{ext}\f$
        summed over all medium components with the specified material type at wavelength
        \f$\lambda\f$ in spatial cell with index \f$m\f$. Because no photon packet is provided,
        default values are used for any relevant incoming photon packet properties. For example,
        the radiation is assumed to be unpolarized. */
    double opacityExt(double lambda, int m, MaterialMix::MaterialType type) const;

    /** This function returns the extinction opacity \f$k^\text{ext}=\sum_h k_h^\text{ext}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$. Because no photon packet is provided, default values are used for any relevant
        incoming photon packet properties. For example, the radiation is assumed to be unpolarized.
        */
    double opacityExt(double lambda, int m) const;

    //=============== High-level photon life cycle ===================

private:
    /** This function returns the absorption opacity \f$k^\text{abs}=\sum_h k_h^\text{abs}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$, where applicable taking into account the properties of the specified incoming
        photon packet (for example, its polarization state). */
    double opacityAbs(double lambda, int m, const PhotonPacket* pp) const;

    /** This function returns the scattering opacity \f$k^\text{sca}=\sum_h k_h^\text{sca}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$, where applicable taking into account the properties of the specified incoming
        photon packet (for example, its polarization state). */
    double opacitySca(double lambda, int m, const PhotonPacket* pp) const;

    /** This function returns the extinction opacity \f$k^\text{ext}=\sum_h k_h^\text{ext}\f$
        summed over all medium components at wavelength \f$\lambda\f$ in spatial cell with index
        \f$m\f$, where applicable taking into account the properties of the specified incoming
        photon packet (for example, its polarization state). */
    double opacityExt(double lambda, int m, const PhotonPacket* pp) const;

public:
    /** This function returns the perceived wavelength of the photon packet at the scattering
        interaction distance, taking into account the bulk velocity and Hubble expansion velocity
        in that cell. */
    double perceivedWavelengthForScattering(const PhotonPacket* pp) const;

    /** This function returns the weighted scattering albedo \f[\frac{\sum_h k_h^\text{sca}}
        {\sum_h k_h^\text{ext}}\f] over all medium components at wavelength \f$\lambda\f$ in the
        spatial cell hosting the specified photon packet's scattering event. The opacities are
        calculated at the wavelength perceived by the medium taking into account the bulk velocity
        and Hubble expansion velocity in that cell and taking into account any relevant properties
        of the specified photon packet such as the polarization state. */
    double albedoForScattering(const PhotonPacket* pp) const;

    /** This function calculates the relative weights of the medium components in a scattering
        event, determined by the scattering opacity \f$k_{m,h}^\text{sca}\f$ of the medium
        component \f$h\f$ in the scattering interaction cell \f$m\f$ obtained from the specified
        photon packet. These opacities are calculated at the specified wavelength (which is assumed
        to be the wavelength perceived by the medium in cell \f$m\f$ taking into account the bulk
        velocity and Hubble expansion velocity in that cell) and taking into account any relevant
        properties of the specified photon packet such as the polarization state.

        The resulting weights are normalized to a total of unity and stored in the target array.
        The array is resized appropriately (i.e. to the number of medium components in the
        simulation). The function returns true if normalized weights have been successfully
        calculated, and false if all of the weights are zero (i.e. the photon packet does not
        scatter in this cell). */
    bool weightsForScattering(ShortArray& wv, double lambda, const PhotonPacket* pp) const;

    /** This function calculates the consolidated peel-off photon luminosity and polarization state
        for all medium components with the specified opacity weights and for the given perceived
        wavelength, geometry, and incoming photon packet. The specified placeholder peel-off photon
        packet is then launched using this information so that it is ready for detection by
        instruments. If there are multiple medium components, the contributions to the luminosity
        (and if polarization is enabled, to the other components of the Stokes vector) are weighted
        by the specified relative opacities of the various medium components.

        This function should be called only when scattering events cannot change the wavelength,
        i.e. the hasScatteringDispersion() function returns false for all material mixes in the
        simulation. In this case, a single consolidated peel-off photon packet can be sent to each
        intrument. */
    void peelOffScattering(const ShortArray& wv, double lambda, Direction bfkobs, Direction bfky, PhotonPacket* pp,
                           PhotonPacket* ppp) const;

    /** This function calculates the peel-off photon luminosity, polarization state, and wavelength
        for the specified medium component with the specified opacity weight and for the given
        perceived wavelength, geometry, and incoming photon packet. The specified placeholder
        peel-off photon packet is then launched using this information so that it is ready for
        detection by instruments.

        This function should be used when scattering events may change the wavelength, i.e. the
        hasScatteringDispersion() function returns true for one or more material mixes in the
        simulation. In this case, a seperate peel-off photon packet for each medium component must
        be sent to each intrument. */
    void peelOffScattering(int h, double w, double lambda, Direction bfkobs, Direction bfky, PhotonPacket* pp,
                           PhotonPacket* ppp) const;

    /** This function simulates a random walk scattering event of a photon packet. Most of the
        properties of the photon packet remain unaltered, including the position and the
        luminosity. The properties that change include the number of scattering events experienced
        by the photon packet, which is increased by one, the propagation direction, which is
        generated randomly, the wavelength, which is properly Doppler-shifted for the bulk velocity
        of the medium, and the polarization state, which may be affected by the scattering process.

        If there is only one medium component, the scattering event is governed by the
        corresponding material mix. If there are several components, the function first randomly
        selects a medium component from the list, where the relative weight of each component
        \f$h\f$ is determined by the scattering opacity \f$k_{m,h}^\text{sca}\f$ of the medium
        component in the scattering interaction cell \f$m\f$ obtained from the specified photon
        packet. These opacities are calculated at the wavelength perceived by the medium in cell
        \f$m\f$ taking into account the bulk velocity and Hubble expansion velocity in that cell,
        and taking into account any relevant properties of the incoming photon packet such as the
        polarization state.

        Performing the actual scattering event is delegated to the material mix corresponding to
        the selected medium component in the interaction cell. Refer to the
        MaterialMix::performScattering() function for more information. */
    void simulateScattering(Random* random, PhotonPacket* pp) const;

    /** This function calculates the cumulative extinction optical depth at the end of each path
        segment along a path through the medium system defined by the initial position and
        direction of the specified PhotonPacket object, and stores the results of the calculation
        into the same PhotonPacket object.

        This function is intended for handling random-walk photon packet paths during a photon life
        cycle that \em does use forced-scattering, requiring the optical depth to be calculated for
        the full path until it reaches the model boundary. Moreover, the function is intended for a
        photon life cycle that does \em not use explicit absorption, so it suffices to consider the
        extinction optical depth (without discriminating between scattering and absorption).

        Because the function is at the heart of the photon life cycle, performance is important.
        Firstly, separating the geometric and optical depth calculations seems to be faster,
        probably due to memory access and caching issues. So the function first determines and
        stores the path segments and then calculates and stores the cumulative optical depth at the
        end of each segment. Secondly, the function implements optimized versions for media with
        spatially constant cross sections.

        With the geometric path information given, the function calculates the extinction optical
        depth for each path segment \f$(\Delta s)_m\f$ as it crosses the spatial cell with index
        \f$m\f$ as \f[ \tau_m^\text{ext} = (\Delta s)_m \sum_h k_{m,h}^\text{ext}, \f] where
        \f$k_{m,h}^\text{ext}\f$ is the extinction opacity corresponding to the \f$h\f$'th medium
        component in the cell with index \f$m\f$ and the sum over \f$h\f$ runs over all medium
        components. The opacities \f$k_{m,h}^\text{ext}\f$ are calculated at the wavelength
        perceived by the medium in cell \f$m\f$ taking into account the bulk velocity and Hubble
        expansion velocity in that cell, and taking into account any relevant properties of the
        incoming photon packet such as the polarization state.

        Using these optical depth values per segment, the function determines the cumulative
        optical depths at the segment exit boundaries and stores them into the specified photon
        packet object. Note that the optical depth at entry of the initial segment is equal to zero
        by definition. */
    void setExtinctionOpticalDepths(PhotonPacket* pp) const;

    /** This function calculates the cumulative scattering and absorption optical depths at the end
        of each path segment along a path through the medium system defined by the initial position
        and direction of the specified PhotonPacket object, and stores the results of the
        calculation into the same PhotonPacket object.

        This function is intended for handling random-walk photon packet paths during a photon life
        cycle that \em does use forced-scattering, requiring the optical depths to be calculated
        for the full path until it reaches the model boundary. Moreover, the function is intended
        for a photon life cycle that \em does use explicit absorption. In this case, the next
        interaction point is determined based on the scattering optical depth (as opposed to the
        total extinction optical depth), so we need to calculate and store the scattering and
        absorption optical depths separately.

        Because it is at the heart of the photon life cycle, performance is important. Firstly,
        separating the geometric and optical depth calculations seems to be faster, probably due to
        memory access and caching issues. So the function first determines and stores the path
        segments and then calculates and stores the cumulative optical depths at the end of each
        segment. Secondly, the function implements optimized versions for media with spatially
        constant cross sections.

        With the geometric path information given, the function calculates the scattering and
        absorption optical depths for each path segment \f$(\Delta s)_m\f$ as it crosses the
        spatial cell with index \f$m\f$ as \f[ \tau_m^\text{sca,abs} = (\Delta s)_m \sum_h
        k_{m,h}^\text{sca,abs}, \f] where \f$k_{m,h}^\text{sca,abs}\f$ represent the scattering and
        absorption opacities corresponding to the \f$h\f$'th medium component in the cell with
        index \f$m\f$ and the sum over \f$h\f$ runs over all medium components. The opacities
        \f$k_{m,h}^\text{sca,abs}\f$ are calculated at the wavelength perceived by the medium in
        cell \f$m\f$ taking into account the bulk velocity and Hubble expansion velocity in that
        cell, and taking into account any relevant properties of the incoming photon packet such as
        the polarization state.

        Using these optical depth values per segment, the function determines the cumulative
        optical depths at the segment exit boundaries and stores them into the specified photon
        packet object. Note that the optical depth at entry of the initial segment is equal to zero
        by definition. */
    void setScatteringAndAbsorptionOpticalDepths(PhotonPacket* pp) const;

    /** This function calculates the cumulative extinction optical depth and distance at the end of
        path segments along a path through the medium system defined by the initial position and
        direction of the specified PhotonPacket object until the specified interaction optical
        depth has been reached. The function then interpolates the extinction optical depth and
        distance at the interaction point, stores these values in the photon packet, and returns
        true. If the specified interaction optical depth is never reached within the path, the
        function returns false.

        This function is intended for handling random-walk photon packet paths during a photon life
        cycle that does \em not use forced-scattering. In that case there is no need to calculate
        the complete path, substantially boosting performance in high-optical depth media.
        Moreover, the function is intended for a photon life cycle that does \em not use explicit
        absorption, so it suffices to consider the extinction optical depth (without discriminating
        between scattering and absorption). Because the function is at the heart of the photon life
        cycle, performance is important. Hence it implements optimized versions for media with
        spatially constant cross sections.

        The optical depth for each traversed path segment is calculated as described for the
        setExtinctionOpticalDepths() function, i.e. at the wavelength perceived by the medium in
        the cell being crossed and taking into account any relevant properties of the incoming
        photon packet. */
    bool setInteractionPointUsingExtinction(PhotonPacket* pp, double tauinteract) const;

    /** This function calculates the cumulative scattering optical depth, the cumulative absorption
        optical depth, and the distance at the end of each of the path segments along a path
        through the medium system defined by the initial position and direction of the specified
        PhotonPacket object. The calculation proceeds until the \em scattering optical depth
        reaches the specified interaction optical depth. The function then interpolates the
        scattering optical depth, the cumulative absorption optical depth, and the distance at the
        interaction point, stores these values in the photon packet, and returns true. If the
        specified interaction optical depth is never reached within the path, the function returns
        false.

        This function is intended for handling random-walk photon packet paths during a photon life
        cycle that does \em not use forced-scattering. In that case there is no need to calculate
        the complete path, substantially boosting performance in high-optical depth media.
        Moreover, the function is intended for a photon life cycle that \em does use explicit
        absorption, requiring the separate treatment of scattering and absorption optical depth.
        Because the function is at the heart of the photon life cycle, performance is important.
        Hence it implements optimized versions for media with spatially constant cross sections.

        The optical depth for each traversed path segment is calculated as described for the
        setScatteringAndAbsorptionOpticalDepths() function, i.e. at the wavelength perceived by the
        medium in the cell being crossed and taking into account any relevant properties of the
        incoming photon packet. */
    bool setInteractionPointUsingScatteringAndAbsorption(PhotonPacket* pp, double tauinteract) const;

    /** This function calculates and returns the extinction optical depth along a path through the
        medium system defined by the initial position and direction of the specified PhotonPacket
        object and up to the specified distance.

        This function is intended for handling peel-off photon packets during the photon life
        cycle. Because it is at the heart of the photon life cycle, performance is important. Hence
        the function implements optimized versions for media with spatially constant cross
        sections. Furthermore, the calculation is limited to the specified distance along the path.
        More precisely, all path segments with an entry boundary at a cumulative distance along the
        path smaller than the specified distance are included in the calculation, and any remaining
        segments are skipped.

        The optical depth for each traversed path segment is calculated as described for the
        setExtinctionOpticalDepths() function, i.e. at the wavelength perceived by the medium in
        the cell being crossed and taking into account any relevant properties of the incoming
        photon packet.

        <b>High optical depth</b>

        Assuming that the extinction cross section is always positive, we know that the observable
        weight of a peel-off photon packet will be numerically zero as soon as the cumulative
        optical depth along its path is higher than \f$ \tau_\mathrm{max} = \ln(L/L_\mathrm{min})
        \f$, where \f$L\f$ is the weight at the peel-off interaction site, and \f$L_\mathrm{min}\f$
        is the smallest representable positive double-precision floating point number. Hence this
        function aborts the calculation and returns positive infinity when this happens.

        If the extinction cross section can be negative, this optimization cannot be applied
        because the cumulative optical depth could decrease again further along the path. */
    double getExtinctionOpticalDepth(const PhotonPacket* pp, double distance) const;

    /** This function returns the extinction optical depth at the specified wavelength along a path
        through the medium system, taking into account only medium components with the specified
        material type. The starting position and the direction of the path are taken from the
        specified SpatialGridPath object. This function is intended for use from probes and hence
        is not performance-sensitive.

        The function determines the segments of the path \f$(\Delta s)_m\f$ as it crosses the cells
        with indices \f$m\f$ in the spatial grid and calculates the optical depth along the path as
        \f[ \tau_\text{path} = \sum_m (\Delta s)_m \sum_h k_{m,h}^\text{ext}, \f] where
        \f$k_{m,h}^\text{ext}\f$ is the extinction opacity corresponding to the \f$h\f$'th medium
        component in the cell with index \f$m\f$ at the specified wavelength \f$\lambda\f$ and the
        sum over \f$h\f$ runs only over the medium components with the specified material type.

        Because no photon packet is available, the calculation ignores kinematic effects as well as
        any other photon packet properties such as polarization. */
    double getExtinctionOpticalDepth(const SpatialGridPath* path, double lambda, MaterialMix::MaterialType type) const;

    //=============== Radiation field ===================

public:
    /** This function initializes all values of the primary and/or secondary radiation field info
        tables to zero. In simulation modes that record the radiation field, the function should be
        called before starting a simulation segment (i.e. before a set of photon packets is
        launched). If the \em primary flag is true, the function clears both the primary table and
        the stable secondary table (if present). The stable secondary table is cleared as well so
        that we can use its contents even if no secondary of photon packet segment has been
        launched yet. If the flag is false, the function clears just the temporary secondary table,
        so that the stable secondary table remains available for calculating secondary emission
        spectra. */
    void clearRadiationField(bool primary);

    /** This function adds the specified value of \f$L\,\Delta s\f$ to the radiation field bin
        corresponding to the spatial cell index \f$m\f$ and the wavelength index\f$\ell\f$. If the
        \em primary flag is true, the value is added to the primary table; otherwise it is added to
        the temporary secondary table.

        The addition happens in a thread-safe way, so that this function can be called from
        multiple parallel threads, even for the same spatial/wavelength bin. If any of the indices
        are out of range, undefined behavior results. */
    void storeRadiationField(bool primary, int m, int ell, double Lds);

    /** This function accumulates the radiation field between multiple processes. In simulation
        modes that record the radiation field, the function should be called in serial code after
        finishing a simulation segment (i.e. after a before set of photon packets has been
        launched) and before querying the radiation field's contents. If the \em primary flag is
        true, the primary table is synchronized; otherwise the temporary secondary table is
        synchronized and its contents is copied into the stable secondary table. */
    void communicateRadiationField(bool primary);

    /** This function returns a pair of values specifying the bolometric luminosity absorbed by
        dust media across the complete domain of the spatial grid, respectively using the partial
        radiation field stored in the primary table and the stable secondary table. The bolometric
        absorbed luminosity in each cell is calculated as described for the dustLuminosity()
        function. */
    std::pair<double, double> totalDustAbsorbedLuminosity() const;

private:
    /** This function returns the sum of the values in both the primary and the stable secondary
        radiation field tables at the specified cell and wavelength indices. If a table is not
        present, the value for that table is assumed to be zero. */
    double radiationField(int m, int ell) const;

public:
    /** This function returns an array with the mean radiation field intensity
        \f$(J_\lambda)_{\ell,m}\f$ in the spatial cell with index \f$m\f$ at each of the wavelength
        bins \f$\ell\f$ defined by the wavelength grid returned by the
        Configuration::radiationFieldWLG() function.

        The mean intensity is calculated using \f[ (J_\lambda)_{\ell,m} = \frac{ (L\Delta
        s)_{\ell,m} }{4\pi\,V_m\,(\Delta \lambda)_\ell} \f] where \f$\ell\f$ is the index of the
        wavelength bin, \f$(\Delta \lambda)_\ell\f$ is the wavelength bin width, \f$m\f$ is the
        spatial cell index, \f$V_m\f$ is the volume of the cell, and \f$(L\Delta s)_{\ell,m}\f$ has
        been accumulated over all photon packets contributing to the bin. The resulting mean
        intensity \f$J_\lambda\f$ is expressed as an amount of energy per unit of time, per unit of
        area, per unit of wavelength, and per unit of solid angle. */
    Array meanIntensity(int m) const;

    //=============== Indicative temperature ===================

public:
    /** This function returns an indicative temperature \f${\bar{T}}_{m,h}\f$ for the material in
        the medium component with index \f$h\f$ in the spatial cell with index \f$m\f$. If the cell
        does not contain any material for the medium component, the function returns zero.

        The interpretation of the indicative temperature depends heavily on the material type.
        Refer to the description for the MaterialMix::indicativeTemperature() function in the
        various MaterialMix subclasses for more information. Note that, in any case, the indicative
        temperature usually does not correspond to a physical temperature.

        The calculation of the indicative temperature may depend on the radiation field embedding
        the material in the specified spatial cell, it may simply reflect the medium state for the
        cell as defined in the input model during setup, or it may somehow depend on both. If the
        radiation field is not being tracked in the simulation, and the indicative temperature for
        the requested material type depends on it, the function will return zero. If the radiation
        field \em is being tracked, this function assumes that the radiation field information has
        been accumulated and the communicateRadiationField() function has been called before
        invoking this function; see the meanIntensity() function. If this is not the case, the
        behavior is undefined. */
    double indicativeTemperature(int m, int h) const;

    /** This function returns an indicative temperature for the material of the specified type in
        the spatial cell with index \f$m\f$. It obtains an indicative temperature for each
        component with a material mix of the specified type, and averages these temperatures over
        the medium components weighed by their relative mass in the cell. In formula form, for a
        spatial cell \f$m\f$ with components \f$h\f$ of the specified type, the indicative dust
        temperature is defined as \f[{\bar{T}}_m = \frac{\sum_h \rho_{m,h}\,{\bar{T}}_{m,h}}
        {\sum_h \rho_{m,h}} \f] where \f${\bar{T}}_{m,h}\f$ is the indicative temperature for each
        component \f$h\f$. If the cell does not contain any material of the requested type, the
        function returns zero.

        The interpretation of the indicative temperature depends heavily on the material type.
        Refer to the description for the MaterialMix::indicativeTemperature() function in the
        various MaterialMix subclasses for more information. Note that, in any case, the indicative
        temperature usually does not correspond to a physical temperature.

        The calculation of the indicative temperature may depend on the radiation field embedding
        the material in the specified spatial cell, it may simply reflect the medium state for the
        cell as defined in the input model during setup, or it may somehow depend on both. If the
        radiation field is not being tracked in the simulation, and the indicative temperature for
        the requested material type depends on it, the function will return zero. If the radiation
        field \em is being tracked, this function assumes that the radiation field information has
        been accumulated and the communicateRadiationField() function has been called before
        invoking this function; see the meanIntensity() function. If this is not the case, the
        behavior is undefined. */
    double indicativeTemperature(int m, MaterialMix::MaterialType type) const;

    /** This function returns an indicative dust temperature for the spatial cell with index
        \f$m\f$, or zero if the simulation does not track the radiation field. Also see the
        indicativeTemperature() function for more information.

        An indicative temperature for each dust component is obtained by solving the energy balance
        equation under LTE (local thermal equilibrium) assumptions for a single representative
        grain for each dust mix. The resulting temperatures are averaged over the dust components
        present in the spatial cell, weighed by relative mass in the cell. If the cell does not
        contain any dust, the function returns zero.

        Note that the indicative dust temperature does not correspond to a physical temperature.
        The LTE assumption is almost certainly unjustified for a relevant portion of the dust
        grains (depending on the embedding radiation field), and even when ignoring this problem,
        averaging temperatures over the dust components and over the various grain material types
        and grain sizes within a particular dust mix has no clear-cut physical justification nor
        interpretation. */
    double indicativeDustTemperature(int m) const;

    /** This function returns an indicative electron temperature in the spatial cell with index
        \f$m\f$. This temperature is obtained by averaging the temperature over the electron medium
        components present in the spatial cell, weighed by relative mass in each component. If no
        medium component specifies an electron temperature, this function returns zero. Also see
        the indicativeTemperature() function for more information. */
    double indicativeElectronTemperature(int m) const;

    /** This function returns an indicative gas temperature in the spatial cell with index \f$m\f$.
        This temperature is obtained by averaging the temperature over the gas medium components
        present in the spatial cell, weighed by relative mass in each component. If no medium
        component specifies a gas temperature, this function returns zero. Also see the
        indicativeTemperature() function for more information. */
    double indicativeGasTemperature(int m) const;

    //====================== Emission ======================

public:
    /** This function returns the bolometric luminosity \f$L^\text{abs}_{\text{bol},m}\f$ that has
        been absorbed by dust media in the spatial cell with index \f$m\f$.

        The bolometric luminosity is calculated using \f[ L^\text{abs}_{\text{bol},m} = \sum_\ell
        (k^\text{abs}_\text{type})_{\ell,m} \,(L\Delta s)_{\ell,m} \f] where \f$\ell\f$ runs over
        the wavelengths in the simulation's radiation field wavelength grid, \f$m\f$ is the spatial
        cell index, \f$(k^\text{abs}_\text{dust})_{\ell,m}\f$ is the absorption opacity of the dust
        in the cell, and \f$(L\Delta s)_{\ell,m}\f$ has been accumulated over all photon packets
        contributing to the bin. */
    double dustLuminosity(int m) const;

    /** This function returns the combined emission spectrum for all dust media in the spatial cell
        with index \f$m\f$. Depending on a simulation-wide configuration option, the spectrum is
        calculated under the assumption of local thermal equilibrium (LTE) or taking into account
        the effect of stochastically heated dust grains (NLTE). The returned spectrum is
        discretized on the wavelength grid returned by the Configuration::dustEmissionWLG()
        function and has arbitrary normalization. The caller is responsible for properly
        normalizing the spectrum based on the value returned by the dustLuminosity() function. */
    Array dustEmissionSpectrum(int m) const;

    /** This function returns the continuum emission spectrum in the spatial cell with index
        \f$m\f$ for the medium component with index \f$h\f$. It is intended for use with gas medium
        components that support secondary continuum emission. When invoked for other medium
        components, the behavior of the function is undefined. The returned spectrum is discretized
        on the wavelength grid returned by the MaterialMix::emissionWavelengthGrid() function of
        the material mix associated with the specified medium component. It contains absolute
        (specific) luminosity values that do not need further normalization. */
    Array continuumEmissionSpectrum(int m, int h) const;

    /** This function returns the line emission spectrum in the spatial cell with index \f$m\f$ for
        the medium component with index \f$h\f$. It is intended for use with gas medium components
        that support secondary line emission. When invoked for other medium components, the
        behavior of the function is undefined. The returned values correspond to each of the lines
        returned by the MaterialMix::lineEmissionCenters() function of the material mix associated
        with the specified medium component. The values represent absolute line luminosities that
        do not need further normalization. */
    Array lineEmissionSpectrum(int m, int h) const;

    //=============== Dynamic medium state ===================

private:
    /** This function invokes the dynamic medium state recipes configured for this simulation to
        update the medium state. The function returns true if all recipes have converged, and false
        otherwise. It assumes that the radiation field has been calculated and that at least one
        dynamic medium state recipe has been configured for the simulation. */
    bool updateDynamicStateRecipes();

    /** This function updates the medium state for any media in the simulation with an associated
        material mix that supports a primary or secondary dynamic medium state, depending on the
        value of the specified flag. The function returns true if all updates have converged, and
        false otherwise. It assumes that the radiation field has been calculated and that at least
        one medium component in the simulation supports a dynamic medium state of the requested
        type. */
    bool updateDynamicStateMedia(bool primary);

public:
    /** This function updates the primary dynamic medium state (PDMS) for all spatial cells and
        medium components based on the currently established radiation field. It invokes any
        dynamic medium state recipes (instances of a DynamicStateRecipe subclass) configured for
        this simulation and updates the medium state for any media in the simulation with an
        associated material mix (instances of a MaterialMix subclass) that supports a PDMS. The
        function returns true if all updates have converged, and false otherwise.

        This function assumes that the radiation field has been calculated. */
    bool updatePrimaryDynamicMediumState();

    /** This function updates the secondary dynamic medium state (SDMS) for all spatial cells and
        medium components based on the currently established radiation field. It updates the medium
        state for any media in the simulation with an associated material mix (instances of a
        MaterialMix subclass) that supports a SDMS. The function returns true if all updates have
        converged, and false otherwise.

        This function assumes that the radiation field has been calculated. */
    bool updateSecondaryDynamicMediumState();

    //======================== Data Members ========================

private:
    Configuration* _config{nullptr};

    // relevant for any simulation mode that includes a medium
    int _numCells{0};  // index m
    int _numMedia{0};  // index h
    bool _mixPerCell{false};
    vector<const MaterialMix*> _mixv;  // material mixes; indexed on h, or on m and h if mixPerCell is true
    MediumState _state;                // state info for each cell and each medium

    // cached info relevant for any simulation mode that includes a medium
    vector<int> _dust_hv;  // a list of indices for media components containing dust
    vector<int> _gas_hv;   // a list of indices for media components containing gas
    vector<int> _elec_hv;  // a list of indices for media components containing electrons
    vector<int> _pdms_hv;  // a list of indices for media components with a primary dynamic medium state
    vector<int> _sdms_hv;  // a list of indices for media components with a secondary dynamic medium state

    // relevant for any simulation mode that stores the radiation field
    WavelengthGrid* _wavelengthGrid{0};  // index ell
    // each radiation field table has an entry for each cell and each wavelength (indexed on m,ell)
    // - the sum of rf1 and rf2 represents the stable radiation field to be used as input for regular calculations
    // - rf2c serves as a target for storing the secondary radiation field so that rf1+rf2 remain available for
    //   calculating secondary emission spectra while already shooting photons through the grid
    Table<2> _rf1;   // radiation field from primary sources
    Table<2> _rf2;   // radiation field from secondary sources (copied from _rf2c at the appropriate time)
    Table<2> _rf2c;  // radiation field currently being accumulated from secondary sources

    // relevant for any simulation mode that includes dust emission
    int _numDustEmissionWavelengths{0};
};

////////////////////////////////////////////////////////////////

#endif

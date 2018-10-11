/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSYSTEM_HPP
#define MEDIUMSYSTEM_HPP

#include "SimulationItem.hpp"
#include "Array.hpp"
#include "Medium.hpp"
#include "MaterialMix.hpp"
#include "SpatialGrid.hpp"
#include "Table.hpp"
class PhotonPacket;
class Random;
class WavelengthGrid;

//////////////////////////////////////////////////////////////////////

/** An instance of the MediumSystem class represents a complete medium system, which is the
    superposition of one or more transfer media. Each individual medium represents a spatial
    density distribution and defines the material properties of the medium at each location. While
    the specific material properties may vary with location, the fundamental material type must be
    the same throughout the spatial domain for each medium.

    In addition to the media input model, the MediumSystem class includes the spatial grid that
    tessellates the spatial domain of the simulation into cells, and manages the medium state for
    each spatial cell in this grid.

    TO DO: add more info on managing the medium state. */
class MediumSystem : public SimulationItem
{
    ITEM_CONCRETE(MediumSystem, SimulationItem, "a medium system")
        ATTRIBUTE_TYPE_ALLOWED_IF(MediumSystem, "!NoMedium")

    PROPERTY_ITEM_LIST(media, Medium, "the transfer media")
        ATTRIBUTE_DEFAULT_VALUE(media, "GeometricMedium")
        ATTRIBUTE_REQUIRED_IF(media, "!NoMedium")

    PROPERTY_ITEM(grid, SpatialGrid, "the spatial grid")
        ATTRIBUTE_DEFAULT_VALUE(grid,
                            "Dimension3:PolicyTreeSpatialGrid;Dimension2:Cylinder2DSpatialGrid;Sphere1DSpatialGrid")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function calculates and stores initial state information for each cell, including the
        cell volume and the number density for each medium as defined by the input model. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

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

    /** This function returns the volume of the spatial cell with index \f$m\f$. */
    double volume(int m) const;

    /** This function returns the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in
        spatial cell with index \f$m\f$. If there are multiple media components, the aggregate bulk
        velocity \f${\boldsymbol{v}}\f$ is determined by averaging the respective bulk velocities
        over the corresponding number densities, \f[{\boldsymbol{v}} = \frac{\sum_h n_h
        {\boldsymbol{v}}_h} {\sum_h n_h}.\f] */
    Vec bulkVelocity(int m);

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

    /** This function returns the number density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double numberDensity(int m, int h) const;

    /** This function returns the mass density of the medium component with index \f$h\f$ in
        spatial cell with index \f$m\f$. */
    double massDensity(int m, int h) const;

    /** This function returns the material mix corresponding to the medium component with index
        \f$h\f$ in spatial cell with index \f$m\f$. */
    const MaterialMix* mix(int m, int h) const;

    /** This function randomly returns a material mix corresponding to one of the medium components
        in spatial cell with index \f$m\f$. The sampling is weighted by the scattering opacity
        \f$k=n_h\sigma_h^\text{sca}\f$ at wavelength \f$\lambda\f$ of each medium component with
        index \f$h\f$ in the spatial cell with index \f$m\f$. */
    const MaterialMix* randomMixForScattering(Random* random, double lambda, int m) const;

    /** This function returns the scattering opacity \f$k=n_h\sigma_h^\text{sca}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. */
    double opacitySca(double lambda, int m, int h) const;

    /** This function returns the scattering opacity \f$k=\sum_h n_h\sigma_h^\text{sca}\f$ summed
        over all medium components at wavelength \f$\lambda\f$ in spatial cell with index \f$m\f$.
        */
    double opacitySca(double lambda, int m) const;

    /** This function returns the extinction opacity \f$k=n_h\sigma_h^\text{ext}\f$ at wavelength
        \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with index
        \f$m\f$. */
    double opacityExt(double lambda, int m, int h) const;

    /** This function returns the extinction opacity \f$k=\sum_h n_h\sigma_h^\text{ext}\f$ summed
        over all medium components at wavelength \f$\lambda\f$ in spatial cell with index \f$m\f$.
        */
    double opacityExt(double lambda, int m) const;

    /** This function returns the extinction opacity \f$k=\sum_h n_h\sigma_h^\text{ext}\f$ summed
        over all medium components with the specified material type at wavelength \f$\lambda\f$ in
        spatial cell with index \f$m\f$. */
    double opacityExt(double lambda, int m, MaterialMix::MaterialType type) const;

    /** This function returns the scattering albedo \f$\sigma_h^\text{sca}/\sigma_h^\text{ext}\f$
        at wavelength \f$\lambda\f$ of the medium component with index \f$h\f$ in spatial cell with
        index \f$m\f$. */
    double albedo(double lambda, int m, int h) const;

    /** This function returns the weighted scattering albedo \f[\frac{\sum_h
        n_h\sigma_h^\text{sca}} {\sum_h n_h\sigma_h^\text{ext}}\f] over all medium components at
        wavelength \f$\lambda\f$ in spatial cell with index \f$m\f$. */
    double albedo(double lambda, int m) const;

    /** This function returns the optical depth at the specified wavelength along a path through
        the medium system, taking into account only medium components with the specified material
        type. The starting position and the direction of the path are taken from the specified
        SpatialGridPath object.

        The function first calls the SpatialGrid::path() function to store the geometrical
        information on the path through the spatial grid into the SpatialGridPath object. It then
        calculates the optical depth along the path as \f[
        \tau_\text{path}(\lambda,{\boldsymbol{r}},{\boldsymbol{k}}) = \sum_m (\Delta s)_m \sum_h
        \varsigma_{\lambda}^{\text{ext}}\, n_m, \f] where \f$\varsigma_{\lambda}^{\text{abs}}\f$ is
        the extinction cross section corresponding to the \f$h\f$'th medium component at wavelength
        \f$\lambda\f$ and \f$n_{m,h}\f$ the number density in the cell with index \f$m\f$
        corresponding to the \f$h\f$'th medium component, and where the sum runs only over the
        medium components with the specified material type. */
    double opticalDepth(SpatialGridPath* path, double lambda, MaterialMix::MaterialType type);

    /** This function returns the extinction factor along a path through the medium system, taking
        into account all medium components. The wavelength, the starting position and the direction
        of the path are taken from the specified PhotonPacket object. The path length is limited to
        the specified distance, if this is smaller than the distance to the edge of the spatial
        grid.

        The function first calls the SpatialGrid::path() function to store the geometrical
        information on the path through the spatial grid into the PhotonPacket object. It then
        calculates the optical depth at the specified distance as as \f[
        \tau_\text{path}(\lambda,{\boldsymbol{r}},{\boldsymbol{k},d}) = \sum_m^{s_m<d} (\Delta s)_m
        \sum_h \varsigma_{\lambda_m,h}^{\text{ext}}\, n_m, \f] where
        \f$\varsigma_{\lambda_m,h}^{\text{abs}}\f$ is the extinction cross section corresponding to
        the \f$h\f$'th medium component at wavelength \f$\lambda_m\f$ and \f$n_{m,h}\f$ the number
        density in the cell with index \f$m\f$ corresponding to the \f$h\f$'th medium component.
        The wavelength \f$\lambda_m\f$ is the wavelength perceived by the medium in cell \f$m\f$
        taking into account the bulk velocity in that cell. The sum over the cells is limited to
        the cells that fall inside the specified distance. Finally, the function calculates the
        exinction factor from the optical depth through \f$\exp(-\tau)\f$. */
    double extinctionFactor(PhotonPacket* pp, double distance);

    /** This function calculates the extinction along a path through the media system defined by
        the specified PhotonPacket object, and stores the results of the calculation into the same
        PhotonPacket object.

        The function first calls the SpatialGrid::path() function to store the geometrical
        information on the path through the spatial grid into the photon packet object, using the
        initial position \f${\boldsymbol{r}}\f$ and direction \f${\boldsymbol{k}}\f$ obtained from
        the photon packet.

        With this information given, the function calculates the optical depth for each path
        segment (or equivalently, for each crossed cell \f$m\f$) as \f[ \tau_m(\lambda_m) = (\Delta
        s)_m \sum_h \varsigma_{\lambda_m,h}^{\text{ext}}\, n_m, \f] where
        \f$\varsigma_{\lambda_m,h}^{\text{abs}}\f$ is the extinction cross section corresponding to
        the \f$h\f$'th medium component at wavelength \f$\lambda_m\f$ and \f$n_{m,h}\f$ the number
        density in the cell with index \f$m\f$ corresponding to the \f$h\f$'th medium component.
        The wavelength \f$\lambda_m\f$ is the wavelength perceived by the medium in cell \f$m\f$
        taking into account the bulk velocity in that cell.

        Subsequently, the function calculates the cumulative optical depth values at the segment
        exit boundaries and stores them into the specified photon packet object. Note that the
        optical depth at entry of the initial segment is equal to zero by definition. */
    void fillOpticalDepthInfo(PhotonPacket* pp);

    /** This function initializes the radiation field memory. In simulation modes that record the
        radiation field, the function should be called before starting a simulation segment (i.e.
        before a set of photon packets is launched). The function properly resizes the radiation
        field data table to include a bin for each spatial cell in the simulation and for each bin
        in the wavelength grid returned by the Configuration::radiationFieldWavelengthGrid()
        function. Finally, the function initializes all values to zero. */
    void clearRadiationField();

    /** This function adds the specified value of \f$L\,\Delta s\f$ to the radiation field bin
        corresponding to the spatial cell index \f$m\f$ and the wavelength index\f$\ell\f$. The
        addition happens in a thread-safe way, so that this function can be called from multiple
        parallel threads, even for the same spatial/wavelength bin. If any of the indices are out
        of range, undefined behavior results. */
    void storeRadiationField(int m, int ell, double Lds);

    /** This function accumulates the radiation field between multiple processes. In simulation
        modes that record the radiation field, the function should be called when finishing a
        simulation segment (i.e. after a before set of photon packets has been launched) and before
        querying the radiation field's contents. */
    void communicateRadiationField();

    /** This function returns an array with the mean radiation field intensity
        \f$(J_\lambda)_{\ell,m}\f$ in the spatial cell with index \f$m\f$ at each of the wavelength
        bins \f$\ell\f$ defined by the wavelength grid returned by the
        Configuration::radiationFieldWavelengthGrid() function.

        This function assumes that a set of photon packets have been launched for a particular
        simulation segment, and that radiation field information has been accumulated during their
        life cycles by calling the storeRadiationField() function. Furthermore, the
        communicateRadiationField() function must have been called before invoking this function.
        If this is not the case, the behavior is undefined.

        The mean intensity is calculated using \f[ (J_\lambda)_{\ell,m} = \frac{ (L\Delta
        s)_{\ell,m} }{4\pi\,V_m\,(\Delta \lambda)_\ell} \f] where \f$\ell\f$ is the index of the
        wavelength bin, \f$(\Delta \lambda)_\ell\f$ is the wavelength bin width, \f$m\f$ is the
        spatial cell index, \f$V_m\f$ is the volume of the cell, and \f$(L\Delta s)_{\ell,m}\f$ has
        been accumulated over all photon packets contributing to the bin. The resulting mean
        intensity \f$J_\lambda\f$ is expressed as an amount of energy per unit of time, per unit of
        area, per unit of wavelength, and per unit of solid angle. */
    Array meanIntensity(int m) const;

    //================== Private Types and Functions ====================

private:
    /** This data structure holds the information maintained per cell. */
    struct State1
    {
        double V;                   // volume
        Vec v;                      // bulk velocity
    };

    /** This data structure holds the information maintained per cell and per medium. */
    struct State2
    {
        double n;                   // the number density
        const MaterialMix* mix;     // pointer to the material mix
    };

    /** This function returns a writable reference to the state data structure for the given cell
        index. */
    State1& state(int m) { return _state1v[m]; }

    /** This function returns a read-only reference to the state data structure for the given cell
        index. */
    const State1& state(int m) const { return _state1v[m]; }

    /** This function returns a writable reference to the state data structure for the given cell
        and medium indices. */
    State2& state(int m, int h) { return _state2vv[m*_numMedia+h]; }

    /** This function returns a read-only reference to the state data structure for the given cell
        and medium indices. */
    const State2& state(int m, int h) const { return _state2vv[m*_numMedia+h]; }

    /** This function communicates the cell states between multiple processes after the states have
        been initialized in parallel (i.e. each process initialized a subset of the states). */
    void communicateStates();

    //======================== Data Members ========================

private:
    // relevant for any simulation mode that includes a medium
    int _numCells{0};           // index m
    int _numMedia{0};           // index h
    vector<State1> _state1v;    // state info for each cell (indexed on m)
    vector<State2> _state2vv;   // state info for each cell and each medium (indexed on m,h)

    // relevant for any simulation mode that stores the radiation field
    WavelengthGrid* _wavelengthGrid{0};  // index ell
    Table<2> _radiationField;   // radiation field info for each cell and each wavelength (indexed on m,ell)
};

////////////////////////////////////////////////////////////////

#endif

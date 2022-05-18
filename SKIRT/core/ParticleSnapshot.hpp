/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARTICLESNAPSHOT_HPP
#define PARTICLESNAPSHOT_HPP

#include "Array.hpp"
#include "Snapshot.hpp"
class SmoothingKernel;

////////////////////////////////////////////////////////////////////

/** A ParticleSnapshot object represents the data in a smoothed particle (SPH) snapshot imported
    from a column text file. The class is based on the Snapshot class; it uses the facilities
    offered there to configure and read the snapshot data, and it implements all functions in the
    general Snapshot public interface. In addition it offers functionality that is specific to this
    snapshot type, such as, for example, the ability to use an arbitrary smoothing kernel for the
    particles in the snapshot.

    If the snapshot configuration requires the ability to determine the density at a given spatial
    position, a lot of effort is made to accelerate the density interpolation over a potentially
    large number of smoothed particles. */
class ParticleSnapshot : public Snapshot
{
    //================= Construction - Destruction =================

public:
    /** The default constructor initializes the snapshot in an invalid state; it is provided here
        so that we don't need to expose the inplemetation of the private Particle class. */
    ParticleSnapshot();

    /** The destructor deletes the smoothed particle grid, if it was constructed. */
    ~ParticleSnapshot();

    //========== Reading ==========

public:
    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and finally closes the file by
        calling the base class Snapshot::readAndClose() function. The function also logs some
        statistical information about the import. If the snapshot configuration requires the
        ability to determine the density at a given spatial position, this function builds a data
        structure that accelerates the density interpolation over a potentially large number of
        smoothed particles. */
    void readAndClose() override;

    //========== Configuration ==========

public:
    /** This function sets the smoothing kernel used for interpolating the smoothed particles. This
        function must be called during configuration. There is no default; failing to set the
        smoothing kernel results in undefined behavior. */
    void setSmoothingKernel(const SmoothingKernel* kernel);

    //=========== Interrogation ==========

public:
    /** This function returns the bounding box surrounding all particles lined up with the
        coordinate axes. The function assumes a kernel with finite support to determine the extent
        of each particle. */
    Box extent() const override;

    /** This function returns the number of particles in the snapshot. */
    int numEntities() const override;

    /** This function returns the effective volume of the particle with index \em m. If no density
        policy has been set or no mass information is being imported, or if the index is out of
        range, the behavior is undefined. */
    double volume(int m) const override;

    /** This function returns the effective mass density associated with the particle with index
        \em m. If no density policy has been set or no mass information is being imported, or if
        the index is out of range, the behavior is undefined. */
    double density(int m) const override;

    /** This function returns the mass density represented by the snapshot at a given point
        \f${\bf{r}}\f$, determined by interpolating (conceptually) over all smoothed particles. If
        the point is outside the domain, the function returns zero. If no density policy has been
        set or no mass information is being imported, the behavior is undefined. */
    double density(Position bfr) const override;

    /** This function returns the total mass represented by the snapshot, in other words the sum of
        the masses of all particles. If no density policy has been set or no mass information is
        being imported, the behavior is undefined. */
    double mass() const override;

    /** This function returns the position for the particle with index \em m. If the index is out
        of range, the behavior is undefined. */
    Position position(int m) const override;

    /** This function returns a random position drawn from the smoothing kernel of the particle
        with index \em m. If the index is out of range, the behavior is undefined. */
    Position generatePosition(int m) const override;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. The function first selects
        a random particle from the discrete probability distribution formed by the respective
        particle masses, and then generates a random position from the smoothing kernel of that
        particle. If no density policy has been set or no mass information is being imported, the
        behavior is undefined. */
    Position generatePosition() const override;

protected:
    /** This function returns a reference to an array containing the imported properties (in column
        order) for the particle with index \f$0\le m \le N_\mathrm{ent}-1\f$. If the index is out
        of range, the behavior is undefined. */
    const Array& properties(int m) const override;

    /** This function returns the index \f$0\le m \le N_\mathrm{ent}-1\f$ of the particle whose
        center is nearest to the specified point \f${\bf{r}}\f$, or -1 if the point is outside the
        domain, if there are no particles in the snapshot, or if the search data structures were
        not created. */
    int nearestEntity(Position bfr) const override;

public:
    /** This function replaces the contents of the specified entity collection by the set of
        particles with a smoothing kernel that overlaps the specified point \f${\bf{r}}\f$. The
        weight corresponding to each particle is set to the particle's smoothing kernel value at
        the given point. If the given point does not overlap any particle, the collection will be
        empty. If the search data structures were not created, invoking this function causes
        undefined behavior. */
    void getEntities(EntityCollection& entities, Position bfr) const override;

    /** This function replaces the contents of the specified entity collection by the set of
        particles with a smoothing kernel that overlaps the specified path with starting point
        \f${\bf{r}}\f$ and direction \f${\bf{k}}\f$. The weight of each particle is given by the
        effective length seen by the path as it crosses the particle's smoothing kernel. If the
        path does not overlap any particle, the collection will be empty. If the search data
        structures were not created, invoking this function causes undefined behavior. */
    void getEntities(EntityCollection& entities, Position bfr, Direction bfk) const override;

    //======================== Data Members ========================

private:
    // private classes
    class Particle;
    class ParticleGrid;

    // data members initialized during configuration
    const SmoothingKernel* _kernel{nullptr};

    // data members initialized when reading the input file
    vector<Array> _propv;  // particle properties as imported

    // data members initialized when reading the input file, but only if a density policy has been set
    vector<Particle> _pv;          // compact particle objects in the same order
    ParticleGrid* _grid{nullptr};  // smart grid for locating smoothed particles
    Array _cumrhov;                // cumulative density distribution for particles
    double _mass{0.};              // total effective mass
};

////////////////////////////////////////////////////////////////////

#endif

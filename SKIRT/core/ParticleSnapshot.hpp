/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARTICLESNAPSHOT_HPP
#define PARTICLESNAPSHOT_HPP

#include "Array.hpp"
#include "SmoothedParticle.hpp"
#include "Snapshot.hpp"
class SmoothedParticleGrid;
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

    /** This function returns the position for the particle with index \em m. If the index is out
        of range, the behavior is undefined. */
    Position position(int m) const override;

    /** This function returns the metallicity of the particle with index \em m. If the metallicity
        is not being imported, or the index is out of range, the behavior is undefined. */
    double metallicity(int m) const override;

    /** This function returns the metallicity of the particle centered nearest to the specified
        point \f${\bf{r}}\f$. If the point is outside the domain, the function returns -1. If the
        metallicity is not being imported, the behavior is undefined. */
    double metallicity(Position bfr) const override;

    /** This function returns the temperature of the particle with index \em m. If the temperature
        is not being imported, or the index is out of range, the behavior is undefined. */
    double temperature(int m) const override;

    /** This function returns the temperature of the particle centered nearest to the specified
        point \f${\bf{r}}\f$. If the point is outside the domain, the function returns -1. If the
        temperature is not being imported, the behavior is undefined. */
    double temperature(Position bfr) const override;

    /** This function returns the velocity of the particle with index \em m. If the velocity is not
        being imported, or the index is out of range, the behavior is undefined. */
    Vec velocity(int m) const override;

    /** This function returns the velocity of the particle centered nearest to the specified point
        \f${\bf{r}}\f$. If the point is outside the domain, the function returns zero velocity. If
        the velocity is not being imported, the behavior is undefined. */
    Vec velocity(Position bfr) const override;

    /** This function returns the velocity dispersion of the entity with index \f$0\le m \le
        N_\mathrm{ent}-1\f$. If the velocity dispersion is not being imported, or the index is out
        of range, the behavior is undefined. */
    double velocityDispersion(int m) const override;

    /** This function returns the velocity dispersion of the entity nearest to (or at) the
        specified point \f${\bf{r}}\f$. If the point is outside the domain, the function returns
        zero dispersion. If the velocity dispersion is not being imported, the behavior is
        undefined. */
    double velocityDispersion(Position bfr) const override;

    /** This function returns the magnetic field vector of the particle with index \em m. If the
        magnetic field is not being imported, or the index is out of range, the behavior is
        undefined. */
    Vec magneticField(int m) const override;

    /** This function returns the magnetic field vector of the particle centered nearest to the
        specified point \f${\bf{r}}\f$. If the point is outside the domain, the function returns a
        zero magnetic field. If the magnetic field is not being imported, the behavior is
        undefined. */
    Vec magneticField(Position bfr) const override;

    /** This function stores the parameters of the particle with index \em m into the given array.
        If parameters are not being imported, or the index is out of range, the behavior is
        undefined. */
    void parameters(int m, Array& params) const override;

    /** This function stores the parameters of the particle centered nearest to the specified point
        \f${\bf{r}}\f$ into the given array. If the point is outside the domain, the function
        returns the appropriate number of zero parameter values. If parameters are not being
        imported, the behavior is undefined. */
    void parameters(Position bfr, Array& params) const override;

    /** This function returns the mass density represented by the snapshot at a given point
        \f${\bf{r}}\f$, determined by interpolating (conceptually) over all smoothed particles. If
        the point is outside the domain, the function returns zero. If no density policy has been
        set or no mass information is being imported, the behavior is undefined. */
    double density(Position bfr) const override;

    /** This function returns the total mass represented by the snapshot, in other words the sum of
        the masses of all particles. If no density policy has been set or no mass information is
        being imported, the behavior is undefined. */
    double mass() const override;

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

    //======================== Data Members ========================

private:
    // data members initialized during configuration
    const SmoothingKernel* _kernel{nullptr};

    // data members initialized when reading the input file
    vector<Array> _propv;  // particle properties as imported

    // data members initialized when reading the input file, but only if a density policy has been set
    vector<SmoothedParticle> _pv;          // compact particle objects in the same order
    SmoothedParticleGrid* _grid{nullptr};  // smart grid for locating smoothed particles
    Array _cumrhov;                        // cumulative density distribution for particles
    double _mass{0.};                      // total effective mass
};

////////////////////////////////////////////////////////////////////

#endif

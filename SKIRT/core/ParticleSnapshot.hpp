/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARTICLESNAPSHOT_HPP
#define PARTICLESNAPSHOT_HPP

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
    position, a lot of effort is made to optimize the density interpolation over a potentially
    large number of smoothed particles. */
class ParticleSnapshot : public Snapshot
{
    //========== Reading ==========

public:
    /** This function creates an input file object corresponding to the specified file and opens it
        for reading, as described for the Snapshot::open() function in the base class. It is
        overridden here to configure the import of the position and smoothing length columns before
        a client starts configuring any optional columns. */
    void open(const SimulationItem* item, string filename, string description) override;

    /** This function reads the snapshot data from the input file, honoring the options set through
        the configuration functions, stores the data for later use, and finally closes the file by
        calling the base class Snapshot::readAndClose() function. The function also logs some
        statistical information about the import. */
    void readAndClose() override;

    //========== Configuration ==========

public:
    /** This function sets the smoothing kernel used for interpolating the smoothed particles.
     TO DO: specify defaults and implement arbitary kernels. */
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

    /** This function returns the velocity of the particle with index \em m. If the velocity is not
        being imported, or the index is out of range, the behavior is undefined. */
    Vec velocity(int m) const override;

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
        with index \em m. If the index is out of range, the behavior is undefined. The first
        argument provides the random generator to be used. */
    Position generatePosition(Random* random, int m) const override;

    /** This function returns a random position within the spatial domain of the snapshot, drawn
        from the mass density distribution represented by the snapshot. The function first selects
        a random particle from the discrete probability distribution formed by the respective
        particle masses, and then generates a random position from the smoothing kernel of that
        particle. If no density policy has been set or no mass information is being imported, the
        behavior is undefined. The first argument provides the random generator to be used. */
    Position generatePosition(Random* random) const override;

    //======================== Data Members ========================

private:
    // data members initialized during configuration
    const SmoothingKernel* _kernel{nullptr};

    // data members initialized when reading the input file
};

////////////////////////////////////////////////////////////////////

#endif

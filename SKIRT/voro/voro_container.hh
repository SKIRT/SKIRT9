// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file voro_container.hh
 * \brief Header file for the Voronoi container class. */

#ifndef VORO_CONTAINER_HH
#define VORO_CONTAINER_HH

#include "voro_cell.hh"

namespace voro {

/** \brief Container class representing a particle system in a three-dimensional
 * rectangular box.
 *
 * This class represents a system of particles in a three-dimensional
 * rectangular box. It contains constants setting the size of the
 * underlying subgrid of blocks that forms the basis of the Voronoi cell
 * computations. It also constructs bound tables that are used in the Voronoi
 * cell computation. */
class container {

    // ---- data members ----

public:
    /** The number of blocks in the x direction. */
    const int nx;
    /** The number of blocks in the y direction. */
    const int ny;
    /** The number of blocks in the z direction. */
    const int nz;
    /** A constant, set to the value of nx multiplied by ny, which
     * is used in the routines that step through blocks in
     * sequence. */
    const int nxy;
    /** A constant, set to the value of nx*ny*nz, which is used in
     * the routines that step through blocks in sequence. */
    const int nxyz;
    /** The size of a computational block in the x direction. */
    const double boxx;
    /** The size of a computational block in the y direction. */
    const double boxy;
    /** The size of a computational block in the z direction. */
    const double boxz;
    /** The inverse box length in the x direction. */
    const double xsp;
    /** The inverse box length in the y direction. */
    const double ysp;
    /** The inverse box length in the z direction. */
    const double zsp;
    /** An array to hold the minimum distances associated with the
     * worklists. This array is initialized during container
     * construction, by the initialize_radii() routine. */
    double *mrad;
    /** The pre-computed block worklists. */
    static const unsigned int wl[wl_seq_length*wl_hgridcu];

    /** The minimum x coordinate of the container. */
    const double ax;
    /** The maximum x coordinate of the container. */
    const double bx;
    /** The minimum y coordinate of the container. */
    const double ay;
    /** The maximum y coordinate of the container. */
    const double by;
    /** The minimum z coordinate of the container. */
    const double az;
    /** The maximum z coordinate of the container. */
    const double bz;

    /** This array holds the numerical IDs of each particle in each
     * computational box. */
    int **id;
    /** A two dimensional array holding particle positions. For the
     * derived container_poly class, this also holds particle
     * radii. */
    double **p;
    /** This array holds the number of particles within each
     * computational box of the container. */
    int *co;
    /** This array holds the maximum amount of particle memory for
     * each computational box of the container. If the number of
     * particles in a particular box ever approaches this limit,
     * more is allocated using the add_particle_memory() function.
     */
    int *mem;
    /** The amount of memory in the array structure for each
     * particle. This is set to 3 when the basic class is
     * initialized, so that the array holds (x,y,z) positions. If
     * the container class is initialized as part of the derived
     * class container_poly, then this is set to 4, to also hold
     * the particle radii. */
    const int ps;

    // ---- methods ----

public:
    /** The constructor sets up the geometry of the container.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
     *                       coordinate directions. */
    container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,int nx_,int ny_,int nz_);

    /** The destructor frees the dynamically allocated memory. */
    ~container();

    /** Put a particle into the correct region of the container.
     * \param[in] n the numerical ID of the inserted particle.
     * \param[in] (x,y,z) the position vector of the inserted particle. */
    void put(int n,double x,double y,double z);

private:
    /** This routine takes a particle position vector, tries to remap it into the
     * primary domain. If successful, it computes the region into which it can be
     * stored and checks that there is enough memory within this region to store
     * it.
     * \param[out] ijk the region index.
     * \param[in,out] (x,y,z) the particle position, remapped into the primary
     *                        domain if necessary.
     * \return True if the particle can be successfully placed into the container,
     * false otherwise. */
    bool put_locate_block(int &ijk,double &x,double &y,double &z);

    /** Takes a particle position vector and computes the region index into which
     * it should be stored. If the particle position is not in the primary domain,
     * the routine bails out.
     * \param[out] ijk the region index.
     * \param[in,out] (x,y,z) the particle position, remapped into the primary
     *                        domain if necessary.
     * \return True if the particle can be successfully placed into the container,
     * false otherwise. */
    bool put_remap(int &ijk,double &x,double &y,double &z);

    /** Increase memory for a particular region.
     * \param[in] i the index of the region to reallocate. */
    void add_particle_memory(int i);

    /** A custom int function that returns consistent stepping
     * for negative numbers, so that (-1.5, -0.5, 0.5, 1.5) maps
     * to (-2,-1,0,1).
     * \param[in] a the number to consider.
     * \return The value of the custom int operation. */
    int step_int(double a) {return a<0?int(a)-1:int(a);}

    /** Computes the minimum distance from a subregion to a given block. If this distance
     * is smaller than the value of minr, then it passes
     * \param[in,out] minr a pointer to the current minimum distance. If the distance
     *                     computed in this function is smaller, then this distance is
     *                     set to the new one.
     * \param[out] (xlo,ylo,zlo) the lower coordinates of the subregion being
     *                           considered.
     * \param[out] (xhi,yhi,zhi) the upper coordinates of the subregion being
     *                           considered.
     * \param[in] (ti,tj,tk) the coordinates of the block. */
    void compute_minimum(double &minr,double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi,int ti,int tj,int tk);

public:
    /** Initializes the Voronoi cell prior to a compute_cell
     * operation for a specific particle being carried out by a
     * voro_compute class. The cell is initialized to fill the
     * entire container.
     * \param[in,out] c a reference to a voronoicell object.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within its block.
     * \param[in] (ci,cj,ck) the coordinates of the block in the
     * 			 container coordinate system.
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] (x,y,z) the position of the particle.
     * \param[out] disp a block displacement used internally by the
     *		    compute_cell routine.
     * \return true. */
    bool initialize_cell(cell &c,int ijk,int q,int ci,int cj,int ck,
            int &i,int &j,int &k,double &x,double &y,double &z,int &disp) {
        double x1,x2,y1,y2,z1,z2,*pp=p[ijk]+ps*q;
        x=*(pp++);y=*(pp++);z=*pp;
        x1=ax-x;x2=bx-x;i=ci;
        y1=ay-y;y2=by-y;j=cj;
        z1=az-z;z2=bz-z;k=ck;
        c.init(x1,x2,y1,y2,z1,z2);
        disp=ijk-i-nx*(j+ny*k);
        return true;
    }

    /** Returns the position of a particle currently being computed
     * relative to the computational block that it is within. It is
     * used to select the optimal worklist entry to use.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] (ci,cj,ck) the block that the particle is within.
     * \param[out] (fx,fy,fz) the position relative to the block.
     */
    void frac_pos(double x,double y,double z,double ci,double cj,double ck,
            double &fx,double &fy,double &fz) {
        fx=x-ax-boxx*ci;
        fy=y-ay-boxy*cj;
        fz=z-az-boxz*ck;
    }

    /** Calculates the index of block in the container structure
     * corresponding to given coordinates.
     * \param[in] (ci,cj,ck) the coordinates of the original block
     * 			 in the current computation, relative
     * 			 to the container coordinate system (not used).
     * \param[in] (ei,ej,ek) the displacement of the current block
     * 			 from the original block.
     * \param[in,out] (qx,qy,qz) the periodic displacement that
     * 			     must be added to the particles
     * 			     within the computed block (not used).
     * \param[in] disp a block displacement used internally by the
     * 		    find_voronoi_cell and compute_cell routines.
     * \return The block index. */
    int region_index(int ci,int cj,int ck,int ei,int ej,int ek,double &qx,double &qy,double &qz,int &disp) {
        (void)ci; (void)cj; (void)ck; (void)qx; (void)qy; (void)qz;
        return disp+ei+nx*(ej+ny*ek);
    }

    /** This is called prior to computing a Voronoi cell for a
     * given particle to initialize any required constants.
     * \param[in] ijk the block that the particle is within.
     * \param[in] s the index of the particle within the block. */
    void r_init(int ijk,int s) {
        (void)ijk; (void)s;
    }

    /** Sets a required constant to be used when carrying out a
     * plane bounds check. */
    void r_prime(double rv) {
        (void)rv;
    }

    /** Carries out a radius bounds check.
     * \param[in] crs the radius squared to be tested.
     * \param[in] mrs the current maximum distance to a Voronoi
     *                vertex multiplied by two.
     * \return True if particles at this radius could not possibly
     * cut the cell, false otherwise. */
    bool r_ctest(double crs,double mrs) {return crs>mrs;}

    /** Scales a plane displacement during a plane bounds check.
     * \param[in] lrs the plane displacement.
     * \return The scaled value. */
    double r_cutoff(double lrs) {return lrs;}

    /** Adds the maximum radius squared to a given value.
     * \param[in] rs the value to consider.
     * \return The value with the radius squared added. */
    double r_max_add(double rs) {return rs;}

    /** Subtracts the radius squared of a particle from a given
     * value.
     * \param[in] rs the value to consider.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return The value with the radius squared subtracted. */
    double r_current_sub(double rs,int ijk,int q) {
        (void)ijk; (void)q;
        return rs;
    }

    /** Scales a plane displacement prior to use in the plane cutting
     * algorithm.
     * \param[in] rs the initial plane displacement.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return The scaled plane displacement. */
    double r_scale(double rs,int ijk,int q) {
        (void)ijk; (void)q;
        return rs;
    }

    /** Scales a plane displacement prior to use in the plane
     * cutting algorithm, and also checks if it could possibly cut
     * the cell.
     * \param[in,out] rs the plane displacement to be scaled.
     * \param[in] mrs the current maximum distance to a Voronoi
     *                vertex multiplied by two.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell could possibly cut the cell, false
     * otherwise. */
    bool r_scale_check(double &rs,double mrs,int ijk,int q) {
        (void)ijk; (void)q;
        return rs<mrs;
    }
};

}

#endif

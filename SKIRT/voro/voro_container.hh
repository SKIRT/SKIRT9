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

/** \brief Class containing data structures common across all particle container classes.
 *
 * This class contains constants and data structures that are common across all
 * particle container classes. It contains constants setting the size of the
 * underlying subgrid of blocks that forms the basis of the Voronoi cell
 * computations. It also constructs bound tables that are used in the Voronoi
 * cell computation, and contains a number of routines that are common across
 * all container classes. */
class voro_base {
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
        bool contains_neighbor(const char* format);
        voro_base(int nx_,int ny_,int nz_,double boxx_,double boxy_,double boxz_);
        ~voro_base() {delete [] mrad;}
    protected:
        /** A custom int function that returns consistent stepping
         * for negative numbers, so that (-1.5, -0.5, 0.5, 1.5) maps
         * to (-2,-1,0,1).
         * \param[in] a the number to consider.
         * \return The value of the custom int operation. */
        inline int step_int(double a) {return a<0?int(a)-1:int(a);}
        /** A custom modulo function that returns consistent stepping
         * for negative numbers. For example, (-2,-1,0,1,2) step_mod 2
         * is (0,1,0,1,0).
         * \param[in] (a,b) the input integers.
         * \return The value of a modulo b, consistent for negative
         * numbers. */
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        /** A custom integer division function that returns consistent
         * stepping for negative numbers. For example, (-2,-1,0,1,2)
         * step_div 2 is (-1,-1,0,0,1).
         * \param[in] (a,b) the input integers.
         * \return The value of a div b, consistent for negative
         * numbers. */
        inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
    private:
        void compute_minimum(double &minr,double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi,int ti,int tj,int tk);
};


/** \brief Class for representing a particle system in a three-dimensional
 * rectangular box.
 *
 * This class represents a system of particles in a three-dimensional
 * rectangular box. Any combination of non-periodic and periodic coordinates
 * can be used in the three coordinate directions. The class is not intended
 * for direct use, but instead forms the base of the container and
 * container_poly classes that add specialized routines for computing the
 * regular and radical Voronoi tessellations respectively. It contains routines
 * that are commonly between these two classes, such as those for
 * placing particles within the internal data structure.
 */
class container_base : public voro_base {
    public:
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
        /** A boolean value that determines if the x coordinate in
         * periodic or not. */
        const bool xperiodic;
        /** A boolean value that determines if the y coordinate in
         * periodic or not. */
        const bool yperiodic;
        /** A boolean value that determines if the z coordinate in
         * periodic or not. */
        const bool zperiodic;
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
        container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
                int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,
                int init_mem,int ps_);
        ~container_base();
        void region_count();
        /** Initializes the Voronoi cell prior to a compute_cell
         * operation for a specific particle being carried out by a
         * voro_compute class. The cell is initialized to fill the
         * entire container. For non-periodic coordinates, this is set
         * by the container boundaries. For periodic coordinates, the
         * space is equally divided in either direction from the
         * particle's initial position.
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
        inline bool initialize_voronoicell(cell &c,int ijk,int q,int ci,int cj,int ck,
                int &i,int &j,int &k,double &x,double &y,double &z,int &disp) {
            double x1,x2,y1,y2,z1,z2,*pp=p[ijk]+ps*q;
            x=*(pp++);y=*(pp++);z=*pp;
            if(xperiodic) {x1=-(x2=0.5*(bx-ax));i=nx;} else {x1=ax-x;x2=bx-x;i=ci;}
            if(yperiodic) {y1=-(y2=0.5*(by-ay));j=ny;} else {y1=ay-y;y2=by-y;j=cj;}
            if(zperiodic) {z1=-(z2=0.5*(bz-az));k=nz;} else {z1=az-z;z2=bz-z;k=ck;}
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
         * 			 to the container coordinate system.
         * \param[in] (ei,ej,ek) the displacement of the current block
         * 			 from the original block.
         * \param[in,out] (qx,qy,qz) the periodic displacement that
         * 			     must be added to the particles
         * 			     within the computed block.
         * \param[in] disp a block displacement used internally by the
         * 		    find_voronoi_cell and compute_cell routines.
         * \return The block index. */
        int region_index(int ci,int cj,int ck,int ei,int ej,int ek,double &qx,double &qy,double &qz,int &disp) {
            if(xperiodic) {if(ci+ei<nx) {ei+=nx;qx=-(bx-ax);} else if(ci+ei>=(nx<<1)) {ei-=nx;qx=bx-ax;} else qx=0;}
            if(yperiodic) {if(cj+ej<ny) {ej+=ny;qy=-(by-ay);} else if(cj+ej>=(ny<<1)) {ej-=ny;qy=by-ay;} else qy=0;}
            if(zperiodic) {if(ck+ek<nz) {ek+=nz;qz=-(bz-az);} else if(ck+ek>=(nz<<1)) {ek-=nz;qz=bz-az;} else qz=0;}
            return disp+ei+nx*(ej+ny*ek);
        }
        /** Sums up the total number of stored particles.
         * \return The number of particles. */
        int total_particles() {
            int tp=*co;
            for(int *cop=co+1;cop<co+nxyz;cop++) tp+=*cop;
            return tp;
        }
    protected:
        void add_particle_memory(int i);
        bool put_locate_block(int &ijk,double &x,double &y,double &z);
        bool put_remap(int &ijk,double &x,double &y,double &z);
        bool remap(int &ai,int &aj,int &ak,int &ci,int &cj,int &ck,double &x,double &y,double &z,int &ijk);
};

/** \brief Class containing all of the routines that are specific to computing
 * the regular Voronoi tessellation.
 *
 * The container and container_periodic classes are derived from this class,
 * and during the Voronoi cell computation, these routines are used to create
 * the regular Voronoi tessellation. */
class radius_mono {
    protected:
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

/** \brief Extension of the container_base class for computing regular Voronoi
 * tessellations.
 *
 * This class is an extension of the container_base class that has routines
 * specifically for computing the regular Voronoi tessellation with no
 * dependence on particle radii. */
class container : public container_base, public radius_mono {
    public:
        container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
                int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem);
        void put(int n,double x,double y,double z);
        friend class compute;
};

}

#endif

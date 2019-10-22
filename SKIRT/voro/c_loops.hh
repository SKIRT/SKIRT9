// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file c_loops.hh
 * \brief Header file for the loop classes. */

#ifndef VOROPP_C_LOOPS_HH
#define VOROPP_C_LOOPS_HH

#include "common.hh"

namespace voro {

/** A type associated with a c_loop_subset class, determining what type of
 * geometrical region to loop over. */
enum c_loop_subset_mode {
    sphere,
    box,
    no_check
};

/** \brief A class for storing ordering information when particles are added to
 * a container.
 *
 * When particles are added to a container class, they are sorted into an
 * internal computational grid of blocks. The particle_order class provides a
 * mechanism for remembering which block particles were sorted into. The
 * put routines in the container class have variants that also take a
 * particle_order class. Each time they are called, they will store the block
 * that the particle was sorted into, plus the position of the particle within
 * the block. The particle_order class can used by the c_loop_order class to
 * specifically loop over the particles that have their information stored
 * within it. */
class particle_order {
    public:
        /** A pointer to the array holding the ordering. */
        int *o;
        /** A pointer to the next position in the ordering array in
         * which to store an entry. */
        int *op;
        /** The current memory allocation for the class, set to the
         * number of entries which can be stored. */
        int size;
        /** The particle_order constructor allocates memory to store the
         * ordering information.
         * \param[in] init_size the initial amount of memory to
         *                      allocate. */
        particle_order(int init_size=init_ordering_size)
            : o(new int[init_size<<1]),op(o),size(init_size) {}
        /** The particle_order destructor frees the dynamically allocated
         * memory used to store the ordering information. */
        ~particle_order() {
            delete [] o;
        }
        /** Adds a record to the order, corresponding to the memory
         * address of where a particle was placed into the container.
         * \param[in] ijk the block into which the particle was placed.
         * \param[in] q the position within the block where the
         * 		particle was placed. */
        inline void add(int ijk,int q) {
            if(op==o+size) add_ordering_memory();
            *(op++)=ijk;*(op++)=q;
        }
    private:
        void add_ordering_memory();
};

/** \brief Base class for looping over particles in a container.
 *
 * This class forms the base of all classes that can loop over a subset of
 * particles in a contaner in some order. When initialized, it stores constants
 * about the corresponding container geometry. It also contains a number of
 * routines for interrogating which particle currently being considered by the
 * loop, which are common between all of the derived classes. */
class c_loop_base {
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
        /** The number of floating point numbers per particle in the
         * associated container data structure. */
        const int ps;
        /** A pointer to the particle position information in the
         * associated container data structure. */
        double **p;
        /** A pointer to the particle ID information in the associated
         * container data structure. */
        int **id;
        /** A pointer to the particle counts in the associated
         * container data structure. */
        int *co;
        /** The current x-index of the block under consideration by the
         * loop. */
        int i;
        /** The current y-index of the block under consideration by the
         * loop. */
        int j;
        /** The current z-index of the block under consideration by the
         * loop. */
        int k;
        /** The current index of the block under consideration by the
         * loop. */
        int ijk;
        /** The index of the particle under consideration within the current
         * block. */
        int q;
        /** The constructor copies several necessary constants from the
         * base container class.
         * \param[in] con the container class to use. */
        template<class c_class>
        c_loop_base(c_class &con) : nx(con.nx), ny(con.ny), nz(con.nz),
                        nxy(con.nxy), nxyz(con.nxyz), ps(con.ps),
                        p(con.p), id(con.id), co(con.co) {}
        /** Returns the position vector of the particle currently being
         * considered by the loop.
         * \param[out] (x,y,z) the position vector of the particle. */
        inline void pos(double &x,double &y,double &z) {
            double *pp=p[ijk]+ps*q;
            x=*(pp++);y=*(pp++);z=*pp;
        }
        /** Returns the ID, position vector, and radius of the particle
         * currently being considered by the loop.
         * \param[out] pid the particle ID.
         * \param[out] (x,y,z) the position vector of the particle.
         * \param[out] r the radius of the particle. If no radius
         * 		 information is available the default radius
         * 		 value is returned. */
        inline void pos(int &pid,double &x,double &y,double &z,double &r) {
            pid=id[ijk][q];
            double *pp=p[ijk]+ps*q;
            x=*(pp++);y=*(pp++);z=*pp;
            r=ps==3?default_radius:*(++pp);
        }
        /** Returns the x position of the particle currently being
         * considered by the loop. */
        inline double x() {return p[ijk][ps*q];}
        /** Returns the y position of the particle currently being
         * considered by the loop. */
        inline double y() {return p[ijk][ps*q+1];}
        /** Returns the z position of the particle currently being
         * considered by the loop. */
        inline double z() {return p[ijk][ps*q+2];}
        /** Returns the ID of the particle currently being considered
         * by the loop. */
        inline int pid() {return id[ijk][q];}
};

/** \brief Class for looping over all of the particles in a container.
 *
 * This is one of the simplest loop classes, that scans the computational
 * blocks in order, and scans all the particles within each block in order. */
class c_loop_all : public c_loop_base {
    public:
        /** The constructor copies several necessary constants from the
         * base container class.
         * \param[in] con the container class to use. */
        template<class c_class>
        c_loop_all(c_class &con) : c_loop_base(con) {}
        /** Sets the class to consider the first particle.
         * \return True if there is any particle to consider, false
         * otherwise. */
        inline bool start() {
            i=j=k=ijk=q=0;
            while(co[ijk]==0) if(!next_block()) return false;
            return true;
        }
        /** Finds the next particle to test.
         * \return True if there is another particle, false if no more
         * particles are available. */
        inline bool inc() {
            q++;
            if(q>=co[ijk]) {
                q=0;
                do {
                    if(!next_block()) return false;
                } while(co[ijk]==0);
            }
            return true;
        }
    private:
        /** Updates the internal variables to find the next
         * computational block with any particles.
         * \return True if another block is found, false if there are
         * no more blocks. */
        inline bool next_block() {
            ijk++;
            i++;
            if(i==nx) {
                i=0;j++;
                if(j==ny) {
                    j=0;k++;
                    if(ijk==nxyz) return false;
                }
            }
            return true;
        }
};

}

#endif

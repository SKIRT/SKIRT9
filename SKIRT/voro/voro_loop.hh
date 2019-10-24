// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file voro_loop.hh
 * \brief Header file for the Voronoi loop class. */

#ifndef VORO_LOOP_HH
#define VORO_LOOP_HH

#include "voro_container.hh"

namespace voro {

/** \brief Class for looping over all of the particles in a container.
 *
 * This class iterates over all of the particles in a container, scanning the
 * computational blocks and all the particles within each block in order. */
class loop {

    // ---- data members ----

private:
    /** The number of blocks in the x direction. */
    const int nx;
    /** The number of blocks in the y direction. */
    const int ny;
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

public:
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

    // ---- methods ----

public:
    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    loop(container &con) : nx(con.nx), ny(con.ny),
                    nxyz(con.nxyz), ps(con.ps),
                    p(con.p), id(con.id), co(con.co) {}
    /** Returns the position vector of the particle currently being
     * considered by the loop.
     * \param[out] (x,y,z) the position vector of the particle. */
    void pos(double &x,double &y,double &z) {
        double *pp=p[ijk]+ps*q;
        x=*(pp++);y=*(pp++);z=*pp;
    }
    /** Returns the x position of the particle currently being
     * considered by the loop. */
    double x() {return p[ijk][ps*q];}
    /** Returns the y position of the particle currently being
     * considered by the loop. */
    double y() {return p[ijk][ps*q+1];}
    /** Returns the z position of the particle currently being
     * considered by the loop. */
    double z() {return p[ijk][ps*q+2];}
    /** Returns the ID of the particle currently being considered
     * by the loop. */
    int pid() {return id[ijk][q];}

    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    bool start() {
        i=j=k=ijk=q=0;
        while(co[ijk]==0) if(!next_block()) return false;
        return true;
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    bool inc() {
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
    bool next_block() {
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

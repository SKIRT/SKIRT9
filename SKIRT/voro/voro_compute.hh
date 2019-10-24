// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file voro_compute.hh
 * \brief Header file for the Voronoi compute class. */

#ifndef VORO_COMPUTE_HH
#define VORO_COMPUTE_HH

#include "voro_loop.hh"

namespace voro {

/** \brief Class for carrying out Voronoi cell computations. */
class compute {

    // ---- data members ----

private:
    /** A reference to the container class on which to carry out */
    container &con;
    /** The size of an internal computational block in the x
     * direction. */
    const double boxx;
    /** The size of an internal computational block in the y
     * direction. */
    const double boxy;
    /** The size of an internal computational block in the z
     * direction. */
    const double boxz;
    /** The inverse box length in the x direction, set to
     * nx/(bx-ax). */
    const double xsp;
    /** The inverse box length in the y direction, set to
     * ny/(by-ay). */
    const double ysp;
    /** The inverse box length in the z direction, set to
     * nz/(bz-az). */
    const double zsp;
    /** The number of boxes in the x direction for the searching mask. */
    const int hx;
    /** The number of boxes in the y direction for the searching mask. */
    const int hy;
    /** The number of boxes in the z direction for the searching mask. */
    const int hz;
    /** A constant, set to the value of hx multiplied by hy, which
     * is used in the routines which step through mask boxes in
     * sequence. */
    const int hxy;
    /** A constant, set to the value of hx*hy*hz, which is used in
     * the routines which step through mask boxes in sequence. */
    const int hxyz;
    /** The number of floating point entries to store for each
     * particle. */
    const int ps;
    /** This array holds the numerical IDs of each particle in each
     * computational box. */
    int **id;
    /** A two dimensional array holding particle positions. For the
     * derived container_poly class, this also holds particle
     * radii. */
    double **p;
    /** An array holding the number of particles within each
     * computational box of the container. */
    int *co;

    /** A constant set to boxx*boxx+boxy*boxy+boxz*boxz, which is
     * frequently used in the computation. */
    const double bxsq;
    /** This sets the current value being used to mark tested blocks
     * in the mask. */
    unsigned int mv;
    /** The current size of the search list. */
    int qu_size;
    /** A pointer to the array of worklists. */
    const unsigned int *wl;
    /** An pointer to the array holding the minimum distances
     * associated with the worklists. */
    double *mrad;
    /** This array is used during the cell computation to determine
     * which blocks have been considered. */
    unsigned int *mask;
    /** An array is used to store the queue of blocks to test
     * during the Voronoi cell computation. */
    int *qu;
    /** A pointer to the end of the queue array, used to determine
     * when the queue is full. */
    int *qu_l;

    // ---- methods ----

public:
    /** The constructor initializes constants from the container class, and
     * sets up the mask and queue used for Voronoi computations.
     * \param[in] con_ a reference to the container class to use.
     * \param[in] (hx_,hy_,hz_) the size of the mask to use. */
    compute(container &con_,int hx_,int hy_,int hz_);

    /** This convenience constructor uses the number of blocks in the specified
     * container as the size of the mask to use. Otherwise it operates just like
     * the other constructor.
     * \param[in] con_ a reference to the container class to use. */
    compute(container &con_) : compute(con_, con_.nx, con_.ny, con_.nz) { }

    /** The class destructor frees the dynamically allocated memory
     * for the mask and queue. */
    ~compute() {
        delete [] qu;
        delete [] mask;
    }

    /** This routine computes a Voronoi cell for a single particle in the
     * container. The algorithm constructs the cell by testing over
     * the neighbors of the particle, working outwards until it reaches those
     * particles which could not possibly intersect the cell. For maximum
     * efficiency, this algorithm is divided into three parts. In the first
     * section, the algorithm tests over the blocks which are in the immediate
     * vicinity of the particle, by making use of one of the precomputed worklists.
     * The code then continues to test blocks on the worklist, but also begins to
     * construct a list of neighboring blocks outside the worklist which may need
     * to be test. In the third section, the routine starts testing these
     * neighboring blocks, evaluating whether or not a particle in them could
     * possibly intersect the cell. For blocks that intersect the cell, it tests
     * the particles in that block, and then adds the block neighbors to the list
     * of potential places to consider.
     * \param[in,out] c a reference to a cell object.
     * \param[in] vl a reference to a loop object indicating the particle
     *               for which to compute the cell.
     * \return False if the Voronoi cell was completely removed during the
     *         computation and has zero volume, true otherwise. */
    bool compute_cell(cell &c,loop &vl) { return compute_cell(c,vl.ijk,vl.q,vl.i,vl.j,vl.k); }

private:
    /** This routine computes a Voronoi cell for a single particle in the
     * container. The algorithm constructs the cell by testing over
     * the neighbors of the particle, working outwards until it reaches those
     * particles which could not possibly intersect the cell. For maximum
     * efficiency, this algorithm is divided into three parts. In the first
     * section, the algorithm tests over the blocks which are in the immediate
     * vicinity of the particle, by making use of one of the precomputed worklists.
     * The code then continues to test blocks on the worklist, but also begins to
     * construct a list of neighboring blocks outside the worklist which may need
     * to be test. In the third section, the routine starts testing these
     * neighboring blocks, evaluating whether or not a particle in them could
     * possibly intersect the cell. For blocks that intersect the cell, it tests
     * the particles in that block, and then adds the block neighbors to the list
     * of potential places to consider.
     * \param[in,out] c a reference to a cell object.
     * \param[in] ijk the index of the block that the test particle is in.
     * \param[in] s the index of the particle within the test block.
     * \param[in] (ci,cj,ck) the coordinates of the block that the test particle is
     *                       in relative to the container data structure.
     * \return False if the Voronoi cell was completely removed during the
     *         computation and has zero volume, true otherwise. */
    bool compute_cell(cell &c,int ijk,int s,int ci,int cj,int ck);

    /** This function checks to see whether a particular block can possibly have
     * any intersection with a Voronoi cell, for the case when the closest point
     * from the cell center to the block is at a corner.
     * \param[in,out] c a reference to a Voronoi cell.
     * \param[in] (xl,yl,zl) the relative coordinates of the corner of the block
     *                       closest to the cell center.
     * \param[in] (xh,yh,zh) the relative coordinates of the corner of the block
     *                       furthest away from the cell center.
     * \return False if the block may intersect, true if does not. */
    bool corner_test(cell &c,double xl,double yl,double zl,double xh,double yh,double zh);

    /** This function checks to see whether a particular block can possibly have
     * any intersection with a Voronoi cell, for the case when the closest point
     * from the cell center to the block is on an edge which points along the x
     * direction.
     * \param[in,out] c a reference to a Voronoi cell.
     * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
     *                    block.
     * \param[in] (yl,zl) the relative y and z coordinates of the corner of the
     *                    block closest to the cell center.
     * \param[in] (yh,zh) the relative y and z coordinates of the corner of the
     *                    block furthest away from the cell center.
     * \return False if the block may intersect, true if does not. */
    bool edge_x_test(cell &c,double x0,double yl,double zl,double x1,double yh,double zh);

    /** This function checks to see whether a particular block can possibly have
     * any intersection with a Voronoi cell, for the case when the closest point
     * from the cell center to the block is on an edge which points along the y
     * direction.
     * \param[in,out] c a reference to a Voronoi cell.
     * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
     *                    block.
     * \param[in] (xl,zl) the relative x and z coordinates of the corner of the
     *                    block closest to the cell center.
     * \param[in] (xh,zh) the relative x and z coordinates of the corner of the
     *                    block furthest away from the cell center.
     * \return False if the block may intersect, true if does not. */
    bool edge_y_test(cell &c,double xl,double y0,double zl,double xh,double y1,double zh);

    /** This function checks to see whether a particular block can possibly have
     * any intersection with a Voronoi cell, for the case when the closest point
     * from the cell center to the block is on an edge which points along the z
     * direction.
     * \param[in,out] c a reference to a Voronoi cell.
     * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the block.
     * \param[in] (xl,yl) the relative x and y coordinates of the corner of the
     *                    block closest to the cell center.
     * \param[in] (xh,yh) the relative x and y coordinates of the corner of the
     *                    block furthest away from the cell center.
     * \return False if the block may intersect, true if does not. */
    bool edge_z_test(cell &c,double xl,double yl,double z0,double xh,double yh,double z1);

    /** This function checks to see whether a particular block can possibly have
     * any intersection with a Voronoi cell, for the case when the closest point
     * from the cell center to the block is on a face aligned with the x direction.
     * \param[in,out] c a reference to a Voronoi cell.
     * \param[in] xl the minimum distance from the cell center to the face.
     * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
     *                    block.
     * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
     *                    block.
     * \return False if the block may intersect, true if does not. */
    bool face_x_test(cell &c,double xl,double y0,double z0,double y1,double z1);

    /** This function checks to see whether a particular block can possibly have
     * any intersection with a Voronoi cell, for the case when the closest point
     * from the cell center to the block is on a face aligned with the y direction.
     * \param[in,out] c a reference to a Voronoi cell.
     * \param[in] yl the minimum distance from the cell center to the face.
     * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
     *                    block.
     * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
     *                    block.
     * \return False if the block may intersect, true if does not. */
    bool face_y_test(cell &c,double x0,double yl,double z0,double x1,double z1);

    /** This function checks to see whether a particular block can possibly have
     * any intersection with a Voronoi cell, for the case when the closest point
     * from the cell center to the block is on a face aligned with the z direction.
     * \param[in,out] c a reference to a Voronoi cell.
     * \param[in] zl the minimum distance from the cell center to the face.
     * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
     *                    block.
     * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
     *                    block.
     * \return False if the block may intersect, true if does not. */
    bool face_z_test(cell &c,double x0,double y0,double zl,double x1,double y1);

    /** This routine checks to see whether a point is within a particular distance
     * of a nearby region. If the point is within the distance of the region, then
     * the routine returns true, and computes the maximum distance from the point
     * to the region. Otherwise, the routine returns false.
     * \param[in] (di,dj,dk) the position of the nearby region to be tested,
     *                       relative to the region that the point is in.
     * \param[in] (fx,fy,fz) the displacement of the point within its region.
     * \param[in] (gxs,gys,gzs) the maximum squared distances from the point to the
     *                          sides of its region.
     * \param[out] crs a reference in which to return the maximum distance to the
     *                 region (only computed if the routine returns false).
     * \param[in] mrs the distance to be tested.
     * \return True if the region is further away than mrs, false if the region in
     *         within mrs. */
    bool compute_min_max_radius(int di,int dj,int dk,double fx,double fy,double fz,double gxs,double gys,double gzs,double& crs,double mrs);

    /** Scans the six orthogonal neighbors of a given block and adds them to the
     * queue if they haven't been considered already. It assumes that the queue
     * will definitely have enough memory to add six entries at the end.
     * \param[in] (ei,ej,ek) the block to consider.
     * \param[in,out] qu_e a pointer to the end of the queue. */
    void add_to_mask(int ei,int ej,int ek,int *&qu_e);

    /** Scans a worklist entry and adds any blocks to the queue
     * \param (q) ?
     * \param (mijk) ?
     * \param[in] (ei,ej,ek) the block to consider.
     * \param[in,out] qu_e a pointer to the end of the queue. */
    void scan_bits_mask_add(unsigned int q,unsigned int *mijk,int ei,int ej,int ek,int *&qu_e);

    /** Adds memory to the queue.
     * \param[in,out] qu_s a reference to the queue start pointer.
     * \param[in,out] qu_e a reference to the queue end pointer. */
    void add_list_memory(int*& qu_s,int*& qu_e);

    /** Resets the mask in cases where the mask counter wraps
     * around. */
    void reset_mask() {
        for(unsigned int *mp(mask);mp<mask+hxyz;mp++) *mp=0;
    }
};

}

#endif

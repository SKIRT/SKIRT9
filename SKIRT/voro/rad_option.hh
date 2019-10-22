// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file rad_option.hh
 * \brief Header file for the classes encapsulating functionality for the
 * regular and radical Voronoi tessellations. */

#ifndef VOROPP_RAD_OPTION_HH
#define VOROPP_RAD_OPTION_HH

#include "common.hh"

namespace voro {

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
        inline void r_init(int ijk,int s) {
            (void)ijk; (void)s;
        }
        /** Sets a required constant to be used when carrying out a
         * plane bounds check. */
        inline void r_prime(double rv) {
            (void)rv;
        }
        /** Carries out a radius bounds check.
         * \param[in] crs the radius squared to be tested.
         * \param[in] mrs the current maximum distance to a Voronoi
         *                vertex multiplied by two.
         * \return True if particles at this radius could not possibly
         * cut the cell, false otherwise. */
        inline bool r_ctest(double crs,double mrs) {return crs>mrs;}
        /** Scales a plane displacement during a plane bounds check.
         * \param[in] lrs the plane displacement.
         * \return The scaled value. */
        inline double r_cutoff(double lrs) {return lrs;}
        /** Adds the maximum radius squared to a given value.
         * \param[in] rs the value to consider.
         * \return The value with the radius squared added. */
        inline double r_max_add(double rs) {return rs;}
        /** Subtracts the radius squared of a particle from a given
         * value.
         * \param[in] rs the value to consider.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The value with the radius squared subtracted. */
        inline double r_current_sub(double rs,int ijk,int q) {
            (void)ijk; (void)q;
            return rs;
        }
        /** Scales a plane displacement prior to use in the plane cutting
         * algorithm.
         * \param[in] rs the initial plane displacement.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The scaled plane displacement. */
        inline double r_scale(double rs,int ijk,int q) {
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
        inline bool r_scale_check(double &rs,double mrs,int ijk,int q) {
            (void)ijk; (void)q;
            return rs<mrs;
        }
};

}
#endif

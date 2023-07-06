// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file common.hh
 * \brief Header file for the small helper functions. */

#ifndef VOROPP_COMMON_HH
#define VOROPP_COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "config.hh"

namespace voro {

void check_duplicate(int n,double x,double y,double z,int id,double *qp);

void voro_fatal_error(const char *p,int status);

}

#endif

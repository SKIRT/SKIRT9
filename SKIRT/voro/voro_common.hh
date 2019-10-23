// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file voro_common.hh
 * \brief Header file with standard library includes, configuration constants, and small helper functions. */

#ifndef VORO_COMMON_HH
#define VORO_COMMON_HH

// Define this macro before including the system headers to silence warnings on Windows
// about the C library functions that could cause buffer overruns, such as sprintf()
#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "voro_config.hh"

namespace voro {

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
inline void voro_fatal_error(const char *p,int status) {
    fprintf(stderr,"voro++: %s\n",p);
    exit(status);
}

}

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file common.cc
 * \brief Implementations of the small helper functions. */

#include "common.hh"

namespace voro {

void check_duplicate(int n,double x,double y,double z,int id,double *qp) {
    double dx=*qp-x,dy=qp[1]-y,dz=qp[2]-z;
    if(dx*dx+dy*dy+dz*dz<1e-10) {
        printf("Duplicate: %d (%g,%g,%g) matches %d (%g,%g,%g)\n",n,x,y,z,id,*qp,qp[1],qp[2]);
        exit(1);
    }
}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void voro_fatal_error(const char *p,int status) {
    fprintf(stderr,"voro++: %s\n",p);
    exit(status);
}

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* safe_fopen(const char *filename,const char *mode) {
    FILE *fp=fopen(filename,mode);
    if(fp==NULL) {
        fprintf(stderr,"voro++: Unable to open file '%s'\n",filename);
        exit(VOROPP_FILE_ERROR);
    }
    return fp;
}

}

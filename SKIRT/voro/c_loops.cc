// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file c_loops.cc
 * \brief Function implementations for the loop classes. */

#include "c_loops.hh"

namespace voro {


/** Extends the memory available for storing the ordering. */
void particle_order::add_ordering_memory() {
    int *no=new int[size<<2],*nop=no,*opp=o;
    while(opp<op) *(nop++)=*(opp++);
    delete [] o;
    size<<=1;o=no;op=nop;
}

}

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file container.cc
 * \brief Function implementations for the container and related classes. */

#include "container.hh"

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *			    coordinate directions.
 * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
 *                                               container is periodic in each
 *                                               coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container_base::container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
        int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem,int ps_)
    : voro_base(nx_,ny_,nz_,(bx_-ax_)/nx_,(by_-ay_)/ny_,(bz_-az_)/nz_),
    ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
    xperiodic(xperiodic_), yperiodic(yperiodic_), zperiodic(zperiodic_),
    id(new int*[nxyz]), p(new double*[nxyz]), co(new int[nxyz]), mem(new int[nxyz]), ps(ps_) {
    int l;
    for(l=0;l<nxyz;l++) co[l]=0;
    for(l=0;l<nxyz;l++) mem[l]=init_mem;
    for(l=0;l<nxyz;l++) id[l]=new int[init_mem];
    for(l=0;l<nxyz;l++) p[l]=new double[ps*init_mem];
}

/** The container destructor frees the dynamically allocated memory. */
container_base::~container_base() {
    int l;
    for(l=0;l<nxyz;l++) delete [] p[l];
    for(l=0;l<nxyz;l++) delete [] id[l];
    delete [] id;
    delete [] p;
    delete [] co;
    delete [] mem;
}

/** The class constructor sets up the geometry of container.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
 *                                               container is periodic in each
 *                                               coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block. */
container::container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
    int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem)
    : container_base(ax_,bx_,ay_,by_,az_,bz_,nx_,ny_,nz_,xperiodic_,yperiodic_,zperiodic_,init_mem,3) {}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container::put(int n,double x,double y,double z) {
    int ijk;
    if(put_locate_block(ijk,x,y,z)) {
        id[ijk][co[ijk]]=n;
        double *pp=p[ijk]+3*co[ijk]++;
        *(pp++)=x;*(pp++)=y;*pp=z;
    }
}

/** This routine takes a particle position vector, tries to remap it into the
 * primary domain. If successful, it computes the region into which it can be
 * stored and checks that there is enough memory within this region to store
 * it.
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base::put_locate_block(int &ijk,double &x,double &y,double &z) {
    if(put_remap(ijk,x,y,z)) {
        if(co[ijk]==mem[ijk]) add_particle_memory(ijk);
        return true;
    }
    return false;
}

/** Takes a particle position vector and computes the region index into which
 * it should be stored. If the container is periodic, then the routine also
 * maps the particle position to ensure it is in the primary domain. If the
 * container is not periodic, the routine bails out.
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base::put_remap(int &ijk,double &x,double &y,double &z) {
    int l;

    ijk=step_int((x-ax)*xsp);
    if(xperiodic) {l=step_mod(ijk,nx);x+=boxx*(l-ijk);ijk=l;}
    else if(ijk<0||ijk>=nx) return false;

    int j=step_int((y-ay)*ysp);
    if(yperiodic) {l=step_mod(j,ny);y+=boxy*(l-j);j=l;}
    else if(j<0||j>=ny) return false;

    int k=step_int((z-az)*zsp);
    if(zperiodic) {l=step_mod(k,nz);z+=boxz*(l-k);k=l;}
    else if(k<0||k>=nz) return false;

    ijk+=nx*j+nxy*k;
    return true;
}

/** Takes a position vector and attempts to remap it into the primary domain.
 * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
 *                       with (0,0,0) corresponding to the primary domain.
 * \param[out] (ci,cj,ck) the index of the block that the position vector is
 *                        within, once it has been remapped.
 * \param[in,out] (x,y,z) the position vector to consider, which is remapped
 *                        into the primary domain during the routine.
 * \param[out] ijk the block index that the vector is within.
 * \return True if the particle is within the container or can be remapped into
 * it, false if it lies outside of the container bounds. */
inline bool container_base::remap(int &ai,int &aj,int &ak,int &ci,int &cj,int &ck,double &x,double &y,double &z,int &ijk) {
    ci=step_int((x-ax)*xsp);
    if(ci<0||ci>=nx) {
        if(xperiodic) {ai=step_div(ci,nx);x-=ai*(bx-ax);ci-=ai*nx;}
        else return false;
    } else ai=0;

    cj=step_int((y-ay)*ysp);
    if(cj<0||cj>=ny) {
        if(yperiodic) {aj=step_div(cj,ny);y-=aj*(by-ay);cj-=aj*ny;}
        else return false;
    } else aj=0;

    ck=step_int((z-az)*zsp);
    if(ck<0||ck>=nz) {
        if(zperiodic) {ak=step_div(ck,nz);z-=ak*(bz-az);ck-=ak*nz;}
        else return false;
    } else ak=0;

    ijk=ci+nx*cj+nxy*ck;
    return true;
}

/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
void container_base::add_particle_memory(int i) {
    int l,nmem=mem[i]<<1;

    // Carry out a check on the memory allocation size
    if(nmem>max_particle_memory)
        voro_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);

    // Allocate new memory and copy in the contents of the old arrays
    int *idp=new int[nmem];
    for(l=0;l<co[i];l++) idp[l]=id[i][l];
    double *pp=new double[ps*nmem];
    for(l=0;l<ps*co[i];l++) pp[l]=p[i][l];

    // Update pointers and delete old arrays
    mem[i]=nmem;
    delete [] id[i];id[i]=idp;
    delete [] p[i];p[i]=pp;
}

}

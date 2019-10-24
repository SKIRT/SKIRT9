// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/* \file voro_cell.hh
 * \brief Header file for the Voronoi cell class. */

#ifndef VORO_CELL_HH
#define VORO_CELL_HH

#include "voro_common.hh"

namespace voro {

/** \brief A class representing a single Voronoi cell with neighbor information.
 *
 * This class represents a single Voronoi cell, as a collection of vertices
 * that are connected by edges. The class contains routines for initializing
 * the Voronoi cell to be simple shapes such as a box, tetrahedron, or octahedron.
 * It the contains routines for recomputing the cell based on cutting it
 * by a plane, which forms the key routine for the Voronoi cell computation.
 * It contains numerous routine for computing statistics about the Voronoi cell.
 */
class cell {

    // ---- data members ----

private:
    /** This holds the current size of the arrays ed and nu, which
     * hold the vertex information. If more vertices are created
     * than can fit in this array, then it is dynamically extended
     * using the add_memory_vertices routine. */
    int current_vertices;
    /** This holds the current maximum allowed order of a vertex,
     * which sets the size of the mem, mep, and mec arrays. If a
     * vertex is created with more vertices than this, the arrays
     * are dynamically extended using the add_memory_vorder routine.
     */
    int current_vertex_order;
    /** This sets the size of the main delete stack. */
    int current_delete_size;
    /** This sets the size of the auxiliary delete stack. */
    int current_delete2_size;
    /** This sets the total number of vertices in the current cell.
     */
    int p;
    /** This is the index of particular point in the cell, which is
     * used to start the tracing routines for plane intersection
     * and cutting. These routines will work starting from any
     * point, but it's often most efficient to start from the last
     * point considered, since in many cases, the cell construction
     * algorithm may consider many planes with similar vectors
     * concurrently. */
    int up;
    /** This is a two dimensional array that holds information
     * about the edge connections of the vertices that make up the
     * cell. The two dimensional array is not allocated in the
     * usual method. To account for the fact the different vertices
     * have different orders, and thus require different amounts of
     * storage, the elements of ed[i] point to one-dimensional
     * arrays in the mep[] array of different sizes.
     *
     * More specifically, if vertex i has order m, then ed[i]
     * points to a one-dimensional array in mep[m] that has 2*m+1
     * entries. The first m elements hold the neighboring edges, so
     * that the jth edge of vertex i is held in ed[i][j]. The next
     * m elements hold a table of relations which is redundant but
     * helps speed up the computation. It satisfies the relation
     * ed[ed[i][j]][ed[i][m+j]]=i. The final entry holds a back
     * pointer, so that ed[i+2*m]=i. The back pointers are used
     * when rearranging the memory. */
    int **ed;
    /** This array holds the order of the vertices in the Voronoi
     * cell. This array is dynamically allocated, with its current
     * size held by current_vertices. */
    int *nu;
    /** This in an array with size 3*current_vertices for holding
     * the positions of the vertices. */
    double *pts;

    /** This a one dimensional array that holds the current sizes
     * of the memory allocations for them mep array.*/
    int *mem;
    /** This is a one dimensional array that holds the current
     * number of vertices of order p that are stored in the mep[p]
     * array. */
    int *mec;
    /** This is a two dimensional array for holding the information
     * about the edges of the Voronoi cell. mep[p] is a
     * one-dimensional array for holding the edge information about
     * all vertices of order p, with each vertex holding 2*p+1
     * integers of information. The total number of vertices held
     * on mep[p] is stored in mem[p]. If the space runs out, the
     * code allocates more using the add_memory() routine. */
    int **mep;

    /** This is the delete stack, used to store the vertices which
     * are going to be deleted during the plane cutting procedure.
     */
    int *ds,*stacke;
    /** This is the auxiliary delete stack, which has size set by
     * current_delete2_size. */
    int *ds2,*stacke2;
    /** This stores the current memory allocation for the marginal
     * cases. */
    int current_marginal;
    /** This stores the total number of marginal points which are
     * currently in the buffer. */
    int n_marg;
    /** This array contains a list of the marginal points, and also
     * the outcomes of the marginal tests. */
    int *marg;
    /** The x coordinate of the normal vector to the test plane. */
    double px;
    /** The y coordinate of the normal vector to the test plane. */
    double py;
    /** The z coordinate of the normal vector to the test plane. */
    double pz;
    /** The magnitude of the normal vector to the test plane. */
    double prsq;

    /** This two dimensional array holds the neighbor information
     * associated with each vertex. mne[p] is a one dimensional
     * array which holds all of the neighbor information for
     * vertices of order p. */
    int **mne;
    /** This is a two dimensional array that holds the neighbor
     * information associated with each vertex. ne[i] points to a
     * one-dimensional array in mne[nu[i]]. ne[i][j] holds the
     * neighbor information associated with the jth edge of vertex
     * i. It is set to the ID number of the plane that made the
     * face that is clockwise from the jth edge. */
    int **ne;

    int *paux1;
    int *paux2;

    // ---- methods ----

public:
    /** Constructs a Voronoi cell and sets up the initial memory. */
    cell();

    /** The destructor deallocates all the dynamic memory. */
    ~cell();

    /** Initializes a Voronoi cell as a rectangular box with the given dimensions.
     * It sets up the edge and vertex information, and then
     * sets up the neighbor information, with initial faces being assigned ID
     * numbers from -1 to -6.
     * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
     * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
     * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
    void init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);

    /** Calculates the volume of the Voronoi cell, by decomposing the cell into
     * tetrahedra extending outward from the zeroth vertex, whose volumes are
     * evaluated using a scalar triple product.
     * \return A floating point number holding the calculated volume. */
    double volume();

    /** Computes the maximum radius squared of a vertex from the center of the
     * cell. It can be used to determine when enough particles have been testing an
     * all planes that could cut the cell have been considered.
     * \return The maximum radius squared of a vertex.*/
    double max_radius_squared();

    /** Calculates the centroid of the Voronoi cell, by decomposing the cell into
     * tetrahedra extending outward from the zeroth vertex.
     * \param[out] (cx,cy,cz) references to floating point numbers in which to
     *                        pass back the centroid vector. */
    void centroid(double &cx,double &cy,double &cz);

    /** Returns a vector of the vertex vectors using the local coordinate system.
     * \param[out] v the vector to store the results in. */
    void vertices(std::vector<double> &v);

    /** Returns a vector of the vertex vectors in the global coordinate system.
     * \param[out] v the vector to store the results in.
     * \param[in] (x,y,z) the position vector of the particle in the global
     *                    coordinate system. */
    void vertices(double x,double y,double z,std::vector<double> &v);

    /** For each face, this routine outputs a bracketed sequence of numbers
     * containing a list of all the vertices that make up that face.
     * \param[out] v the vector to store the results in. */
    void face_vertices(std::vector<int> &v);

    /** Returns a list of IDs of neighboring particles
     * corresponding to each face.
     * \param[out] v a reference to a vector in which to return the
     *               results. If no neighbor information is
     *               available, a blank vector is returned. */
    void neighbors(std::vector<int> &v);

    /** Cuts the Voronoi cell by a particle whose center is at a
     * separation of (x,y,z) from the cell center. The value of rsq
     * should be initially set to \f$x^2+y^2+z^2\f$.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \param[in] p_id the plane ID (for neighbor tracking only).
     * \return False if the plane cut deleted the cell entirely,
     * true otherwise. */
    bool nplane(double x,double y,double z,double rsq,int p_id);

    /** This routine tests to see whether the cell intersects a plane by starting
     * from the guess point up. If up intersects, then it immediately returns true.
     * Otherwise, it calls the plane_intersects_track() routine.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \return False if the plane does not intersect the plane, true if it does. */
    bool plane_intersects(double x,double y,double z,double rsq);

    /** This routine tests to see if a cell intersects a plane. It first tests a
     * random sample of approximately sqrt(p)/4 points. If any of those are
     * intersect, then it immediately returns true. Otherwise, it takes the closest
     * point and passes that to plane_intersect_track() routine.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \return False if the plane does not intersect the plane, true if it does. */
    bool plane_intersects_guess(double x,double y,double z,double rsq);

private:
    /** This routine tests to see if a cell intersects a plane, by tracing over the cell from
     * vertex to vertex, starting at up. It is meant to be called either by plane_intersects()
     * or plane_intersects_track(), when those routines cannot immediately resolve the case.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \param[in] g the distance of up from the plane.
     * \return False if the plane does not intersect the plane, true if it does. */
    bool plane_intersects_track(double x,double y,double z,double rsq,double g);

    /** This is a simple function for picking out the index
     * of the next edge counterclockwise at the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return 0 if a=nu[p]-1, or a+1 otherwise. */
    int cycle_up(int a,int p) {return a==nu[p]-1?0:a+1;}

    /** This is a simple function for picking out the index
     * of the next edge clockwise from the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return nu[p]-1 if a=0, or a-1 otherwise. */
    int cycle_down(int a,int p) {return a==0?nu[p]-1:a-1;}

    /** Several routines in the class that gather cell-based statistics internally
     * track their progress by flipping edges to negative so that they know what
     * parts of the cell have already been tested. This function resets them back
     * to positive. When it is called, it assumes that every edge in the routine
     * should have already been flipped to negative, and it bails out with an
     * internal error if it encounters a positive edge. */
    void reset_edges();

    /** Increases the memory storage for a particular vertex order, by increasing
     * the size of the of the corresponding mep array. If the arrays already exist,
     * their size is doubled; if they don't exist, then new ones of size
     * init_n_vertices are allocated. The routine also ensures that the pointers in
     * the ed array are updated, by making use of the back pointers. For the cases
     * where the back pointer has been temporarily overwritten in the marginal
     * vertex code, the auxiliary delete stack is scanned to find out how to update
     * the ed value. If neighbor
     * tracking is turned on, then the routine also reallocates the corresponding mne
     * array.
     * \param[in] i the order of the vertex memory to be increased.
     * \param (stackp2) ?  */
    void add_memory(int i,int *stackp2);

    /** Doubles the maximum number of vertices allowed, by reallocating the ed, nu,
     * and pts arrays. If the allocation exceeds the absolute maximum set in
     * max_vertices, then the routine exits with a fatal error. If neighbor tracking
     * is turned on, then the routine
     * also reallocates the ne array. */
    void add_memory_vertices();

    /** Doubles the maximum allowed vertex order, by reallocating mem, mep, and mec
     * arrays. If the allocation exceeds the absolute maximum set in
     * max_vertex_order, then the routine causes a fatal error. If neighbor tracking
     * is turned on, then the routine also reallocates the mne array. */
    void add_memory_vorder();

    /** Doubles the size allocation of the main delete stack. If the allocation
     * exceeds the absolute maximum set in max_delete_size, then routine causes a
     * fatal error. */
    void add_memory_ds(int *&stackp);

    /** Doubles the size allocation of the auxiliary delete stack. If the
     * allocation exceeds the absolute maximum set in max_delete2_size, then the
     * routine causes a fatal error. */
    void add_memory_ds2(int *&stackp2);

    /** Order one vertices can potentially be created during the order two collapse
     * routine. This routine keeps removing them until there are none left.
     * \return False if the vertex removal was unsuccessful, indicative of the cell
     *         having zero volume and disappearing; true if the vertex removal was
     *         successful. */
    bool collapse_order1();

    /** During the creation of a new facet in the plane routine, it is possible
     * that some order two vertices may arise. This routine removes them.
     * Suppose an order two vertex joins c and d. If there's a edge between
     * c and d already, then the order two vertex is just removed; otherwise,
     * the order two vertex is removed and c and d are joined together directly.
     * It is possible this process will create order two or order one vertices,
     * and the routine is continually run until all of them are removed.
     * \return False if the vertex removal was unsuccessful, indicative of the cell
     *         reducing to zero volume and disappearing; true if the vertex removal
     *         was successful. */
    bool collapse_order2();

    /** This routine deletes the kth edge of vertex j and reorganizes the memory.
     * If the neighbor computation is enabled, we also have to supply an handedness
     * flag to decide whether to preserve the plane on the left or right of the
     * connection.
     * \return False if a zero order vertex was formed, indicative of the cell
     *         disappearing; true if the vertex removal was successful. */
    bool delete_connection(int j,int k,bool hand);

    /** Starting from a point within the current cutting plane, this routine attempts
     * to find an edge to a point outside the cutting plane. This prevents the plane
     * routine from .
     * \param[in,out] up */
    bool search_for_outside_edge(int &up);

    /** Adds a point to the auxiliary delete stack if it is not already there.
     * \param[in] lp the index of the point to add.
     * \param[in,out] stackp2 a pointer to the end of the stack entries. */
    void add_to_stack(int lp,int *&stackp2);

    /** Checks to see if a given vertex is inside, outside or within the test
     * plane. If the point is far away from the test plane, the routine immediately
     * returns whether it is inside or outside. If the routine is close the the
     * plane and within the specified tolerance, then the special check_marginal()
     * routine is called.
     * \param[in] n the vertex to test.
     * \param[out] ans the result of the scalar product used in evaluating the
     *                 location of the point.
     * \return -1 if the point is inside the plane, 1 if the point is outside the
     *         plane, or 0 if the point is within the plane. */
    int m_test(int n,double &ans);

    /** Checks to see if a given vertex is inside, outside or within the test
     * plane, for the case when the point has been detected to be very close to the
     * plane. The routine ensures that the returned results are always consistent
     * with previous tests, by keeping a table of any marginal results. The routine
     * first sees if the vertex is in the table, and if it finds a previously
     * computed result it uses that. Otherwise, it computes a result for this
     * vertex and adds it the table.
     * \param[in] n the vertex to test.
     * \param[in] ans the result of the scalar product used in evaluating
     *                the location of the point.
     * \return -1 if the point is inside the plane, 1 if the point is outside the
     *         plane, or 0 if the point is within the plane. */
    int check_marginal(int n,double &ans);

    void n_allocate(int i,int m) {mne[i]=new int[m*i];}
    void n_add_memory_vertices(int i) {
        int **pp=new int*[i];
        for(int j=0;j<current_vertices;j++) pp[j]=ne[j];
        delete [] ne;ne=pp;
    }
    void n_add_memory_vorder(int i) {
        int **p2=new int*[i];
        for(int j=0;j<current_vertex_order;j++) p2[j]=mne[j];
        delete [] mne;mne=p2;
    }
    void n_set_pointer(int p,int n) {
        ne[p]=mne[n]+n*mec[n];
    }
    void n_copy(int a,int b,int c,int d) {ne[a][b]=ne[c][d];}
    void n_set(int a,int b,int c) {ne[a][b]=c;}
    void n_set_aux1(int k) {paux1=mne[k]+k*mec[k];}
    void n_copy_aux1(int a,int b) {paux1[b]=ne[a][b];}
    void n_copy_aux1_shift(int a,int b) {paux1[b]=ne[a][b+1];}
    void n_set_aux2_copy(int a,int b) {
        paux2=mne[b]+b*mec[b];
        for(int i=0;i<b;i++) ne[a][i]=paux2[i];
    }
    void n_copy_pointer(int a,int b) {ne[a]=ne[b];}
    void n_set_to_aux1(int j) {ne[j]=paux1;}
    void n_set_to_aux2(int j) {ne[j]=paux2;}
    void n_allocate_aux1(int i) {paux1=new int[i*mem[i]];}
    void n_switch_to_aux1(int i) {delete [] mne[i];mne[i]=paux1;}
    void n_copy_to_aux1(int i,int m) {paux1[m]=mne[i][m];}
    void n_set_to_aux1_offset(int k,int m) {ne[k]=paux1+m;}
};

}

#endif

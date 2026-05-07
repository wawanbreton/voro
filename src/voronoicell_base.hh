// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_VORONOICELL_BASE_HH
#define VOROPP_VORONOICELL_BASE_HH

#include <vector>
#include "common.hh"
#include "config.hh"
namespace voro
{

/** \brief A class representing a single Voronoi cell.
 *
 * This class represents a single Voronoi cell, as a collection of vertices
 * that are connected by edges. The class contains routines for initializing
 * the Voronoi cell to be simple shapes such as a box, tetrahedron, or octahedron.
 * It the contains routines for recomputing the cell based on cutting it
 * by a plane, which forms the key routine for the Voronoi cell computation.
 * It contains numerous routine for computing statistics about the Voronoi cell,
 * and it can output the cell in several formats.
 *
 * This class is not intended for direct use, but forms the base of the
 * voronoicell and voronoicell_neighbor classes, which extend it based on
 * whether neighboring particle ID information needs to be tracked. */
class voronoicell_base
{
public:
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
    /** This sets the size of the extra search stack. */
    int current_xsearch_size;
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
    int** ed;
    /** This array holds the order of the vertices in the Voronoi
     * cell. This array is dynamically allocated, with its current
     * size held by current_vertices. */
    int* nu;
    unsigned int* mask;
    /** This in an array with size 3*current_vertices for holding
     * the positions of the vertices. */
    double* pts;
    double tol;
    double tol_cu;
    double big_tol;
    voronoicell_base(double max_len_sq);
    ~voronoicell_base();
    void init_base(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    void init_octahedron_base(double l);
    void init_tetrahedron_base(
        double x0,
        double y0,
        double z0,
        double x1,
        double y1,
        double z1,
        double x2,
        double y2,
        double z2,
        double x3,
        double y3,
        double z3);
    void translate(double x, double y, double z);
    void draw_pov(double x, double y, double z, FILE* fp = stdout);
    /** Outputs the cell in POV-Ray format, using cylinders for edges
     * and spheres for vertices, to a given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    inline void draw_pov(double x, double y, double z, const char* filename)
    {
        FILE* fp = safe_fopen(filename, "w");
        draw_pov(x, y, z, fp);
        fclose(fp);
    };
    void draw_pov_mesh(double x, double y, double z, FILE* fp = stdout);
    /** Outputs the cell in POV-Ray format as a mesh2 object to a
     * given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    inline void draw_pov_mesh(double x, double y, double z, const char* filename)
    {
        FILE* fp = safe_fopen(filename, "w");
        draw_pov_mesh(x, y, z, fp);
        fclose(fp);
    }
    void draw_gnuplot(double x, double y, double z, FILE* fp = stdout);
    /** Outputs the cell in Gnuplot format a given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    inline void draw_gnuplot(double x, double y, double z, const char* filename)
    {
        FILE* fp = safe_fopen(filename, "w");
        draw_gnuplot(x, y, z, fp);
        fclose(fp);
    }
    double volume();
    double max_radius_squared();
    double total_edge_distance();
    double surface_area();
    void centroid(double& cx, double& cy, double& cz);
    int number_of_faces();
    int number_of_edges();
    void vertex_orders(std::vector<int>& v);
    void output_vertex_orders(FILE* fp = stdout);
    void vertices(std::vector<double>& v);
    void output_vertices(FILE* fp = stdout);
    void vertices(double x, double y, double z, std::vector<double>& v);
    void output_vertices(double x, double y, double z, FILE* fp = stdout);
    void solid_angles(std::vector<double>& v);
    void face_areas(std::vector<double>& v);
    void minkowski(double r, double& ar, double& vo);
    /** Outputs the solid angles of the faces.
     * \param[in] fp the file handle to write to. */
    inline void output_solid_angles(FILE* fp = stdout)
    {
        std::vector<double> v;
        solid_angles(v);
        voro_print_vector(v, fp);
    }
    /** Outputs the areas of the faces.
     * \param[in] fp the file handle to write to. */
    inline void output_face_areas(FILE* fp = stdout)
    {
        std::vector<double> v;
        face_areas(v);
        voro_print_vector(v, fp);
    }
    void face_orders(std::vector<int>& v);
    /** Outputs a list of the number of sides of each face.
     * \param[in] fp the file handle to write to. */
    inline void output_face_orders(FILE* fp = stdout)
    {
        std::vector<int> v;
        face_orders(v);
        voro_print_vector(v, fp);
    }
    void face_freq_table(std::vector<int>& v);
    /** Outputs a */
    inline void output_face_freq_table(FILE* fp = stdout)
    {
        std::vector<int> v;
        face_freq_table(v);
        voro_print_vector(v, fp);
    }
    void face_vertices(std::vector<int>& v);
    /** Outputs the */
    inline void output_face_vertices(FILE* fp = stdout)
    {
        std::vector<int> v;
        face_vertices(v);
        voro_print_face_vertices(v, fp);
    }
    void face_perimeters(std::vector<double>& v);
    /** Outputs a list of the perimeters of each face.
     * \param[in] fp the file handle to write to. */
    inline void output_face_perimeters(FILE* fp = stdout)
    {
        std::vector<double> v;
        face_perimeters(v);
        voro_print_vector(v, fp);
    }
    void normals(std::vector<double>& v);
    /** Outputs a list of the perimeters of each face.
     * \param[in] fp the file handle to write to. */
    inline void output_normals(FILE* fp = stdout)
    {
        std::vector<double> v;
        normals(v);
        voro_print_positions(v, fp);
    }
    /** Outputs a custom string of information about the Voronoi
     * cell to a file. It assumes the cell is at (0,0,0) and has a
     * the default_radius associated with it.
     * \param[in] format the custom format string to use.
     * \param[in] fp the file handle to write to. */
    inline void output_custom(const char* format, FILE* fp = stdout)
    {
        output_custom(format, 0, 0, 0, 0, default_radius, fp);
    }
    void output_custom(const char* format, int i, double x, double y, double z, double r, FILE* fp = stdout);
    template<class vc_class>
    bool nplane(vc_class& vc, double x, double y, double z, double rsq, int p_id);
    bool plane_intersects(double x, double y, double z, double rsq);
    bool plane_intersects_guess(double x, double y, double z, double rsq);
    void construct_relations();
    void check_relations();
    void check_duplicates();
    void print_edges();
    /** Returns a list of IDs of neighboring particles
     * corresponding to each face.
     * \param[out] v a reference to a vector in which to return the
     *               results. If no neighbor information is
     *               available, a blank vector is returned. */
    virtual void neighbors(std::vector<int>& v)
    {
        v.clear();
    }
    /** This is a virtual function that is overridden by a routine
     * to print a list of IDs of neighboring particles
     * corresponding to each face. By default, when no neighbor
     * information is available, the routine does nothing.
     * \param[in] fp the file handle to write to. */
    virtual void output_neighbors(FILE* fp = stdout)
    {
    }
    /** This a virtual function that is overridden by a routine to
     * print the neighboring particle IDs for a given vertex. By
     * default, when no neighbor information is available, the
     * routine does nothing.
     * \param[in] i the vertex to consider. */
    virtual void print_edges_neighbors(int i) {};
    /** This is a simple inline function for picking out the index
     * of the next edge counterclockwise at the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return 0 if a=nu[p]-1, or a+1 otherwise. */
    inline int cycle_up(int a, int p)
    {
        return a == nu[p] - 1 ? 0 : a + 1;
    }
    /** This is a simple inline function for picking out the index
     * of the next edge clockwise from the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return nu[p]-1 if a=0, or a-1 otherwise. */
    inline int cycle_down(int a, int p)
    {
        return a == 0 ? nu[p] - 1 : a - 1;
    }

protected:
    /** This a one dimensional array that holds the current sizes
     * of the memory allocations for them mep array.*/
    int* mem;
    /** This is a one dimensional array that holds the current
     * number of vertices of order p that are stored in the mep[p]
     * array. */
    int* mec;
    /** This is a two dimensional array for holding the information
     * about the edges of the Voronoi cell. mep[p] is a
     * one-dimensional array for holding the edge information about
     * all vertices of order p, with each vertex holding 2*p+1
     * integers of information. The total number of vertices held
     * on mep[p] is stored in mem[p]. If the space runs out, the
     * code allocates more using the add_memory() routine. */
    int** mep;
    inline void reset_edges();
    template<class vc_class>
    void check_memory_for_copy(vc_class& vc, voronoicell_base* vb);
    void copy(voronoicell_base* vb);

private:
    /** This is the delete stack, used to store the vertices which
     * are going to be deleted during the plane cutting procedure.
     */
    int *ds, *stackp, *stacke;
    /** This is the auxiliary delete stack, which has size set by
     * current_delete2_size. */
    int *ds2, *stackp2, *stacke2;
    /** This is the extra search stack. */
    int *xse, *stackp3, *stacke3;
    unsigned int maskc;
    /** The x coordinate of the normal vector to the test plane. */
    double px;
    /** The y coordinate of the normal vector to the test plane. */
    double py;
    /** The z coordinate of the normal vector to the test plane. */
    double pz;
    /** The magnitude of the normal vector to the test plane. */
    double prsq;
    template<class vc_class>
    void add_memory(vc_class& vc, int i);
    template<class vc_class>
    void add_memory_vertices(vc_class& vc);
    template<class vc_class>
    void add_memory_vorder(vc_class& vc);
    void add_memory_ds();
    void add_memory_ds2();
    void add_memory_xse();
    bool failsafe_find(int& lp, int& ls, int& us, double& l, double& u);
    template<class vc_class>
    bool create_facet(vc_class& vc, int lp, int ls, double l, int us, double u, int p_id);
    template<class vc_class>
    bool collapse_order1(vc_class& vc);
    template<class vc_class>
    inline bool collapse_order2(vc_class& vc);
    template<class vc_class>
    bool delete_connection(vc_class& vc, int j, int k, bool hand);
    inline bool search_for_outside_edge(int& up);
    inline void add_to_stack(int sc2, int lp);
    inline void reset_mask()
    {
        for (int i = 0; i < current_vertices; i++)
            mask[i] = 0;
        maskc = 4;
    }
    inline bool search_downward(unsigned int& uw, int& lp, int& ls, int& us, double& l, double& u);
    bool definite_max(int& lp, int& ls, double& l, double& u, unsigned int& uw);
    inline bool search_upward(unsigned int& lw, int& lp, int& ls, int& us, double& l, double& u);
    bool definite_min(int& lp, int& us, double& l, double& u, unsigned int& lw);
    inline void minkowski_contrib(int i, int k, int m, double r, double& ar, double& vo);
    void minkowski_edge(double x0, double r1, double s1, double r2, double s2, double r, double& ar, double& vo);
    void minkowski_formula(double x0, double y0, double z0, double r, double& ar, double& vo);
    inline bool plane_intersects_track(double x, double y, double z, double rs, double g);
    inline void normals_search(std::vector<double>& v, int i, int j, int k);
    inline bool search_edge(int l, int& m, int& k);
    inline unsigned int m_test(int n, double& ans);
    inline unsigned int m_testx(int n, double& ans);
    unsigned int m_calc(int n, double& ans);
    inline void flip(int tp)
    {
        ed[tp][nu[tp] << 1] = -1 - ed[tp][nu[tp] << 1];
    }
    int check_marginal(int n, double& ans);
    friend class voronoicell;
    friend class voronoicell_neighbor;
};

} // namespace voro

#endif

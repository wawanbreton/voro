// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_VORONOICELL_NEIGHBOR_HH
#define VOROPP_VORONOICELL_NEIGHBOR_HH

#include "voronoicell.hh"
namespace voro
{

/** \brief Extension of the voronoicell_base class to represent a Voronoi cell
 * with neighbor information.
 *
 * This class is an extension of the voronoicell_base class, in cases when the
 * IDs of neighboring particles associated with each face of the Voronoi cell.
 * It contains additional data structures mne and ne for storing this
 * information. */
class voronoicell_neighbor : public voronoicell_base
{
public:
    using voronoicell_base::nplane;
    /** This two dimensional array holds the neighbor information
     * associated with each vertex. mne[p] is a one dimensional
     * array which holds all of the neighbor information for
     * vertices of order p. */
    int** mne;
    /** This is a two dimensional array that holds the neighbor
     * information associated with each vertex. ne[i] points to a
     * one-dimensional array in mne[nu[i]]. ne[i][j] holds the
     * neighbor information associated with the jth edge of vertex
     * i. It is set to the ID number of the plane that made the
     * face that is clockwise from the jth edge. */
    int** ne;
    voronoicell_neighbor()
        : voronoicell_base(default_length * default_length)
    {
        memory_setup();
    }
    voronoicell_neighbor(double max_len_sq_)
        : voronoicell_base(max_len_sq_)
    {
        memory_setup();
    }
    template<class c_class>
    voronoicell_neighbor(c_class& con)
        : voronoicell_base(con.max_len_sq)
    {
        memory_setup();
    }
    ~voronoicell_neighbor();
    void operator=(voronoicell& c);
    void operator=(voronoicell_neighbor& c);
    /** Cuts the Voronoi cell by a particle whose center is at a
     * separation of (x,y,z) from the cell center. The value of rsq
     * should be initially set to \f$x^2+y^2+z^2\f$.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \param[in] p_id the plane ID (for neighbor tracking only).
     * \return False if the plane cut deleted the cell entirely,
     * true otherwise. */
    inline bool nplane(double x, double y, double z, double rsq, int p_id)
    {
        return nplane(*this, x, y, z, rsq, p_id);
    }
    /** This routine calculates the modulus squared of the vector
     * before passing it to the main nplane() routine with full
     * arguments.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \param[in] p_id the plane ID (for neighbor tracking only).
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool nplane(double x, double y, double z, int p_id)
    {
        double rsq = x * x + y * y + z * z;
        return nplane(*this, x, y, z, rsq, p_id);
    }
    /** This version of the plane routine just makes up the plane
     * ID to be zero. It will only be referenced if neighbor
     * tracking is enabled.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \param[in] rsq the modulus squared of the vector.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z, double rsq)
    {
        return nplane(*this, x, y, z, rsq, 0);
    }
    /** Cuts a Voronoi cell using the influence of a particle at
     * (x,y,z), first calculating the modulus squared of this
     * vector before passing it to the main nplane() routine. Zero
     * is supplied as the plane ID, which will be ignored unless
     * neighbor tracking is enabled.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z)
    {
        double rsq = x * x + y * y + z * z;
        return nplane(*this, x, y, z, rsq, 0);
    }
    void init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    void init_octahedron(double l);
    void init_tetrahedron(
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
    void check_facets();
    virtual void neighbors(std::vector<int>& v);
    virtual void print_edges_neighbors(int i);
    virtual void output_neighbors(FILE* fp = stdout)
    {
        std::vector<int> v;
        neighbors(v);
        voro_print_vector(v, fp);
    }

private:
    int* paux1;
    int* paux2;
    void memory_setup();
    inline void n_allocate(int i, int m)
    {
        mne[i] = new int[m * i];
    }
    inline void n_add_memory_vertices(int i)
    {
        int** pp = new int*[i];
        for (int j = 0; j < current_vertices; j++)
            pp[j] = ne[j];
        delete[] ne;
        ne = pp;
    }
    inline void n_add_memory_vorder(int i)
    {
        int** p2 = new int*[i];
        for (int j = 0; j < current_vertex_order; j++)
            p2[j] = mne[j];
        delete[] mne;
        mne = p2;
    }
    inline void n_set_pointer(int p, int n)
    {
        ne[p] = mne[n] + n * mec[n];
    }
    inline void n_copy(int a, int b, int c, int d)
    {
        ne[a][b] = ne[c][d];
    }
    inline void n_set(int a, int b, int c)
    {
        ne[a][b] = c;
    }
    inline void n_set_aux1(int k)
    {
        paux1 = mne[k] + k * mec[k];
    }
    inline void n_copy_aux1(int a, int b)
    {
        paux1[b] = ne[a][b];
    }
    inline void n_copy_aux1_shift(int a, int b)
    {
        paux1[b] = ne[a][b + 1];
    }
    inline void n_set_aux2_copy(int a, int b)
    {
        paux2 = mne[b] + b * mec[b];
        for (int i = 0; i < b; i++)
            ne[a][i] = paux2[i];
    }
    inline void n_copy_pointer(int a, int b)
    {
        ne[a] = ne[b];
    }
    inline void n_set_to_aux1(int j)
    {
        ne[j] = paux1;
    }
    inline void n_set_to_aux2(int j)
    {
        ne[j] = paux2;
    }
    inline void n_allocate_aux1(int i)
    {
        paux1 = new int[i * mem[i]];
    }
    inline void n_switch_to_aux1(int i)
    {
        delete[] mne[i];
        mne[i] = paux1;
    }
    inline void n_copy_to_aux1(int i, int m)
    {
        paux1[m] = mne[i][m];
    }
    inline void n_set_to_aux1_offset(int k, int m)
    {
        ne[k] = paux1 + m;
    }
    friend class voronoicell_base;
};

} // namespace voro

#endif

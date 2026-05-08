// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_VORONOICELL_HH
#define VOROPP_VORONOICELL_HH

#include "voronoicell_base.hh"
namespace voro
{

/** \brief Extension of the voronoicell_base class to represent a Voronoi
 * cell without neighbor information.
 *
 * This class is an extension of the voronoicell_base class, in cases when
 * is not necessary to track the IDs of neighboring particles associated
 * with each face of the Voronoi cell. */
class voronoicell : public voronoicell_base
{
public:
    using voronoicell_base::nplane;
    voronoicell()
        : voronoicell_base(default_length * default_length)
    {
    }
    voronoicell(double max_len_sq_)
        : voronoicell_base(max_len_sq_)
    {
    }
    template<class c_class>
    voronoicell(c_class& con)
        : voronoicell_base(con.max_len_sq)
    {
    }
    /** Copies the information from another voronoicell class into
     * this class, extending memory allocation if necessary.
     * \param[in] c the class to copy. */
    inline void operator=(voronoicell& c)
    {
        voronoicell_base* vb((voronoicell_base*)&c);
        check_memory_for_copy(*this, vb);
        copy(vb);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] rsq the modulus squared of the vector.
     * \param[in] p_id the plane ID, ignored for this case where no
     *                 neighbor tracking is enabled.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool nplane(double x, double y, double z, double rsq, int p_id)
    {
        return nplane(*this, x, y, z, rsq, 0);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] p_id the plane ID, ignored for this case where no
     *                 neighbor tracking is enabled.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool nplane(double x, double y, double z, int p_id)
    {
        double rsq = x * x + y * y + z * z;
        return nplane(*this, x, y, z, rsq, 0);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] rsq the modulus squared of the vector.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z, double rsq)
    {
        return nplane(*this, x, y, z, rsq, 0);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z)
    {
        double rsq = x * x + y * y + z * z;
        return nplane(*this, x, y, z, rsq, 0);
    }
    /** Initializes the Voronoi cell to be rectangular box with the
     * given dimensions.
     * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
     * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
     * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
    inline void init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
    {
        init_base(xmin, xmax, ymin, ymax, zmin, zmax);
    }
    /** Initializes the cell to be an octahedron with vertices at
     * (l,0,0), (-l,0,0), (0,l,0), (0,-l,0), (0,0,l), and (0,0,-l).
     * \param[in] l a parameter setting the size of the octahedron.
     */
    inline void init_octahedron(double l)
    {
        init_octahedron_base(l);
    }
    /** Initializes the cell to be a tetrahedron.
     * \param[in] (x0,y0,z0) the coordinates of the first vertex.
     * \param[in] (x1,y1,z1) the coordinates of the second vertex.
     * \param[in] (x2,y2,z2) the coordinates of the third vertex.
     * \param[in] (x3,y3,z3) the coordinates of the fourth vertex.
     */
    inline void init_tetrahedron(
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
        double z3)
    {
        init_tetrahedron_base(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    }
    void init_l_shape();

private:
    inline void n_allocate(int i, int m) {};
    inline void n_add_memory_vertices(int i) {};
    inline void n_add_memory_vorder(int i) {};
    inline void n_set_pointer(int p, int n) {};
    inline void n_copy(int a, int b, int c, int d) {};
    inline void n_set(int a, int b, int c) {};
    inline void n_set_aux1(int k) {};
    inline void n_copy_aux1(int a, int b) {};
    inline void n_copy_aux1_shift(int a, int b) {};
    inline void n_set_aux2_copy(int a, int b) {};
    inline void n_copy_pointer(int a, int b) {};
    inline void n_set_to_aux1(int j) {};
    inline void n_set_to_aux2(int j) {};
    inline void n_allocate_aux1(int i) {};
    inline void n_switch_to_aux1(int i) {};
    inline void n_copy_to_aux1(int i, int m) {};
    inline void n_set_to_aux1_offset(int k, int m) {};
    inline void n_neighbors(std::vector<int>& v)
    {
        v.clear();
    };
    friend class voronoicell_base;
};

} // namespace voro

#endif

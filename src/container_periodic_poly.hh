// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_CONTAINER_PERIODIC_POLY_HH
#define VOROPP_CONTAINER_PERIODIC_POLY_HH

#include "container_periodic_base.hh"
#include "particle_order.hh"
#include "c_loop_all_periodic.hh"
#include "radius_poly.hh"
#include "v_compute.hh"
#include "voronoicell.hh"
#include "voronoicell_neighbor.hh"
namespace voro
{

/** \brief Extension of the container_periodic_base class for computing radical
 * Voronoi tessellations.
 *
 * This class is an extension of container_periodic_base that has routines
 * specifically for computing the radical Voronoi tessellation that depends
 * on the particle radii. */
class container_periodic_poly : public container_periodic_base, public radius_poly
{
public:
    container_periodic_poly(
        double bx_,
        double bxy_,
        double by_,
        double bxz_,
        double byz_,
        double bz_,
        int nx_,
        int ny_,
        int nz_,
        int init_mem_);
    void clear();
    void put(int n, double x, double y, double z, double r);
    void put(int n, double x, double y, double z, double r, int& ai, int& aj, int& ak);
    void put(particle_order& vo, int n, double x, double y, double z, double r);
    void import(FILE* fp = stdin);
    void import(particle_order& vo, FILE* fp = stdin);
    /** Imports a list of particles from an open file stream into
     * the container_poly class. Entries of five numbers (Particle
     * ID, x position, y position, z position, radius) are searched
     * for. If the file cannot be successfully read, then the
     * routine causes a fatal error.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(const char* filename)
    {
        FILE* fp = safe_fopen(filename, "r");
        import(fp);
        fclose(fp);
    }
    /** Imports a list of particles from an open file stream into
     * the container_poly class. Entries of five numbers (Particle
     * ID, x position, y position, z position, radius) are searched
     * for. In addition, the order in which particles are read is
     * saved into an ordering class. If the file cannot be
     * successfully read, then the routine causes a fatal error.
     * \param[in,out] vo the ordering class to use.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(particle_order& vo, const char* filename)
    {
        FILE* fp = safe_fopen(filename, "r");
        import(vo, fp);
        fclose(fp);
    }
    void compute_all_cells();
    double sum_cell_volumes();
    /** Dumps particle IDs, positions and radii to a file.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template<class c_loop>
    void draw_particles(c_loop& vl, FILE* fp)
    {
        double* pp;
        if (vl.start())
            do
            {
                pp = p[vl.ijk] + 4 * vl.q;
                fprintf(fp, "%d %g %g %g %g\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2], pp[3]);
            } while (vl.inc());
    }
    /** Dumps all of the particle IDs, positions and radii to a
     * file.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles(FILE* fp = stdout)
    {
        c_loop_all_periodic vl(*this);
        draw_particles(vl, fp);
    }
    /** Dumps all of the particle IDs, positions and radii to a
     * file.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles(const char* filename)
    {
        FILE* fp = safe_fopen(filename, "w");
        draw_particles(fp);
        fclose(fp);
    }
    /** Dumps particle positions in POV-Ray format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template<class c_loop>
    void draw_particles_pov(c_loop& vl, FILE* fp)
    {
        double* pp;
        if (vl.start())
            do
            {
                pp = p[vl.ijk] + 4 * vl.q;
                fprintf(fp, "// id %d\nsphere{<%g,%g,%g>,%g}\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2], pp[3]);
            } while (vl.inc());
    }
    /** Dumps all the particle positions in POV-Ray format.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles_pov(FILE* fp = stdout)
    {
        c_loop_all_periodic vl(*this);
        draw_particles_pov(vl, fp);
    }
    /** Dumps all the particle positions in POV-Ray format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles_pov(const char* filename)
    {
        FILE* fp(safe_fopen(filename, "w"));
        draw_particles_pov(fp);
        fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template<class c_loop>
    void draw_cells_gnuplot(c_loop& vl, FILE* fp)
    {
        voronoicell c(*this);
        double* pp;
        if (vl.start())
            do
                if (compute_cell(c, vl))
                {
                    pp = p[vl.ijk] + ps * vl.q;
                    c.draw_gnuplot(*pp, pp[1], pp[2], fp);
                }
            while (vl.inc());
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_gnuplot(FILE* fp = stdout)
    {
        c_loop_all_periodic vl(*this);
        draw_cells_gnuplot(vl, fp);
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_gnuplot(const char* filename)
    {
        FILE* fp(safe_fopen(filename, "w"));
        draw_cells_gnuplot(fp);
        fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template<class c_loop>
    void draw_cells_pov(c_loop& vl, FILE* fp)
    {
        voronoicell c(*this);
        double* pp;
        if (vl.start())
            do
                if (compute_cell(c, vl))
                {
                    fprintf(fp, "// cell %d\n", id[vl.ijk][vl.q]);
                    pp = p[vl.ijk] + ps * vl.q;
                    c.draw_pov(*pp, pp[1], pp[2], fp);
                }
            while (vl.inc());
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_pov(FILE* fp = stdout)
    {
        c_loop_all_periodic vl(*this);
        draw_cells_pov(vl, fp);
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_pov(const char* filename)
    {
        FILE* fp(safe_fopen(filename, "w"));
        draw_cells_pov(fp);
        fclose(fp);
    }
    /** Computes the Voronoi cells and saves customized information
     * about them.
     * \param[in] vl the loop class to use.
     * \param[in] format the custom output string to use.
     * \param[in] fp a file handle to write to. */
    template<class c_loop>
    void print_custom(c_loop& vl, const char* format, FILE* fp)
    {
        int ijk, q;
        double* pp;
        if (contains_neighbor(format))
        {
            voronoicell_neighbor c(*this);
            if (vl.start())
                do
                    if (compute_cell(c, vl))
                    {
                        ijk = vl.ijk;
                        q = vl.q;
                        pp = p[ijk] + ps * q;
                        c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], pp[3], fp);
                    }
                while (vl.inc());
        }
        else
        {
            voronoicell c(*this);
            if (vl.start())
                do
                    if (compute_cell(c, vl))
                    {
                        ijk = vl.ijk;
                        q = vl.q;
                        pp = p[ijk] + ps * q;
                        c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], pp[3], fp);
                    }
                while (vl.inc());
        }
    }
    /** Computes the Voronoi cell for a particle currently being
     * referenced by a loop class.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] vl the loop class to use.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    template<class v_cell, class c_loop>
    inline bool compute_cell(v_cell& c, c_loop& vl)
    {
        return vc.compute_cell(c, vl.ijk, vl.q, vl.i, vl.j, vl.k);
    }
    /** Computes the Voronoi cell for given particle.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    template<class v_cell>
    inline bool compute_cell(v_cell& c, int ijk, int q)
    {
        int k(ijk / (nx * oy)), ijkt(ijk - (nx * oy) * k), j(ijkt / nx), i(ijkt - j * nx);
        return vc.compute_cell(c, ijk, q, i, j, k);
    }
    /** Computes the Voronoi cell for a ghost particle at a given
     * location.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] (x,y,z) the location of the ghost particle.
     * \param[in] r the radius of the ghost particle.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template<class v_cell>
    inline bool compute_ghost_cell(v_cell& c, double x, double y, double z, double r)
    {
        int ijk;
        put_locate_block(ijk, x, y, z);
        double *pp = p[ijk] + 4 * co[ijk]++, tm = max_radius;
        *(pp++) = x;
        *(pp++) = y;
        *(pp++) = z;
        *pp = r;
        if (r > max_radius)
            max_radius = r;
        bool q = compute_cell(c, ijk, co[ijk] - 1);
        co[ijk]--;
        max_radius = tm;
        return q;
    }
    void print_custom(const char* format, FILE* fp = stdout);
    void print_custom(const char* format, const char* filename);
    bool find_voronoi_cell(double x, double y, double z, double& rx, double& ry, double& rz, int& pid);

private:
    voro_compute<container_periodic_poly> vc;
    friend class voro_compute<container_periodic_poly>;
};

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_CONTAINER_PERIODIC_BASE_HH
#define VOROPP_CONTAINER_PERIODIC_BASE_HH

#include <cstdio>
#include <vector>
#include "common.hh"
#include "config.hh"
#include "unitcell.hh"
#include "v_base.hh"
namespace voro
{

/** \brief Class for representing a particle system in a 3D periodic
 * non-orthogonal periodic domain.
 *
 * This class represents a particle system in a three-dimensional
 * non-orthogonal periodic domain. The domain is defined by three periodicity
 * vectors (bx,0,0), (bxy,by,0), and (bxz,byz,bz) that represent a
 * parallelepiped. Internally, the class stores particles in the box 0<x<bx,
 * 0<y<by, 0<z<bz, and constructs periodic images of particles that displaced
 * by the three periodicity vectors when they are necessary for the
 * computation. The internal memory structure for this class is significantly
 * different from the container_base class in order to handle the dynamic
 * construction of these periodic images.
 *
 * The class is derived from the unitcell class, which encapsulates information
 * about the domain geometry, and the voro_base class, which encapsulates
 * information about the underlying computational grid. */
class container_periodic_base : public unitcell, public voro_base
{
public:
    const double max_len_sq;
    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    int ey;
    /** The lower z index (inclusive) of the primary domain within
     * the block structure. */
    int ez;
    /** The upper y index (exclusive) of the primary domain within
     * the block structure. */
    int wy;
    /** The upper z index (exclusive) of the primary domain within
     * the block structure. */
    int wz;
    /** The total size of the block structure (including images) in
     * the y direction. */
    int oy;
    /** The total size of the block structure (including images) in
     * the z direction. */
    int oz;
    /** The total number of blocks. */
    int oxyz;
    /** This array holds the numerical IDs of each particle in each
     * computational box. */
    int** id;
    /** A two dimensional array holding particle positions. For the
     * derived container_poly class, this also holds particle
     * radii. */
    double** p;
    /** This array holds the number of particles within each
     * computational box of the container. */
    int* co;
    /** This array holds the maximum amount of particle memory for
     * each computational box of the container. If the number of
     * particles in a particular box ever approaches this limit,
     * more is allocated using the add_particle_memory() function.
     */
    int* mem;
    /** An array holding information about periodic image
     * construction at a given location. */
    char* img;
    /** The initial amount of memory to allocate for particles
     * for each block. */
    const int init_mem;
    /** The amount of memory in the array structure for each
     * particle. This is set to 3 when the basic class is
     * initialized, so that the array holds (x,y,z) positions. If
     * the container class is initialized as part of the derived
     * class container_poly, then this is set to 4, to also hold
     * the particle radii. */
    const int ps;
    container_periodic_base(
        double bx_,
        double bxy_,
        double by_,
        double bxz_,
        double byz_,
        double bz_,
        int nx_,
        int ny_,
        int nz_,
        int init_mem_,
        int ps);
    ~container_periodic_base();
    /** Prints all particles in the container, including those that
     * have been constructed in image blocks. */
    inline void print_all_particles()
    {
        int ijk, q;
        for (ijk = 0; ijk < oxyz; ijk++)
            for (q = 0; q < co[ijk]; q++)
                printf("%d %g %g %g\n", id[ijk][q], p[ijk][ps * q], p[ijk][ps * q + 1], p[ijk][ps * q + 2]);
    }
    void region_count();
    /** Initializes the Voronoi cell prior to a compute_cell
     * operation for a specific particle being carried out by a
     * voro_compute class. The cell is initialized to be the
     * pre-computed unit Voronoi cell based on planes formed by
     * periodic images of the particle.
     * \param[in,out] c a reference to a voronoicell object.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within its block.
     * \param[in] (ci,cj,ck) the coordinates of the block in the
     * 			 container coordinate system.
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] (x,y,z) the position of the particle.
     * \param[out] disp a block displacement used internally by the
     *		    compute_cell routine.
     * \return False if the plane cuts applied by walls completely
     * removed the cell, true otherwise. */
    template<class v_cell>
    inline bool initialize_voronoicell(
        v_cell& c,
        int ijk,
        int q,
        int ci,
        int cj,
        int ck,
        int& i,
        int& j,
        int& k,
        double& x,
        double& y,
        double& z,
        int& disp)
    {
        c = unit_voro;
        double* pp = p[ijk] + ps * q;
        x = *(pp++);
        y = *(pp++);
        z = *pp;
        i = nx;
        j = ey;
        k = ez;
        return true;
    }
    /** Initializes parameters for a find_voronoi_cell call within
     * the voro_compute template.
     * \param[in] (ci,cj,ck) the coordinates of the test block in
     * 			 the container coordinate system.
     * \param[in] ijk the index of the test block
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] disp a block displacement used internally by the
     *		    find_voronoi_cell routine (but not needed
     *		    in this instance.) */
    inline void initialize_search(int ci, int cj, int ck, int ijk, int& i, int& j, int& k, int& disp)
    {
        i = nx;
        j = ey;
        k = ez;
    }
    /** Returns the position of a particle currently being computed
     * relative to the computational block that it is within. It is
     * used to select the optimal worklist entry to use.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] (ci,cj,ck) the block that the particle is within.
     * \param[out] (fx,fy,fz) the position relative to the block.
     */
    inline void frac_pos(double x, double y, double z, double ci, double cj, double ck, double& fx, double& fy, double& fz)
    {
        fx = x - boxx * ci;
        fy = y - boxy * (cj - ey);
        fz = z - boxz * (ck - ez);
    }
    /** Calculates the index of block in the container structure
     * corresponding to given coordinates.
     * \param[in] (ci,cj,ck) the coordinates of the original block
     * 			 in the current computation, relative
     * 			 to the container coordinate system.
     * \param[in] (ei,ej,ek) the displacement of the current block
     * 			 from the original block.
     * \param[in,out] (qx,qy,qz) the periodic displacement that
     * 			     must be added to the particles
     * 			     within the computed block.
     * \param[in] disp a block displacement used internally by the
     * 		    find_voronoi_cell and compute_cell routines
     * 		    (but not needed in this instance.)
     * \return The block index. */
    inline int region_index(int ci, int cj, int ck, int ei, int ej, int ek, double& qx, double& qy, double& qz, int& disp)
    {
        int qi = ci + (ei - nx), qj = cj + (ej - ey), qk = ck + (ek - ez);
        int iv(step_div(qi, nx));
        if (iv != 0)
        {
            qx = iv * bx;
            qi -= nx * iv;
        }
        else
            qx = 0;
        create_periodic_image(qi, qj, qk);
        return qi + nx * (qj + oy * qk);
    }
    void create_all_images();
    void check_compartmentalized();

protected:
    void add_particle_memory(int i);
    void put_locate_block(int& ijk, double& x, double& y, double& z);
    void put_locate_block(int& ijk, double& x, double& y, double& z, int& ai, int& aj, int& ak);
    /** Creates particles within an image block by copying them
     * from the primary domain and shifting them. If the given
     * block is aligned with the primary domain in the z-direction,
     * the routine calls the simpler create_side_image routine
     * where the image block may comprise of particles from up to
     * two primary blocks. Otherwise is calls the more complex
     * create_vertical_image where the image block may comprise of
     * particles from up to four primary blocks.
     * \param[in] (di,dj,dk) the coordinates of the image block to
     *                       create. */
    inline void create_periodic_image(int di, int dj, int dk)
    {
        if (di < 0 || di >= nx || dj < 0 || dj >= oy || dk < 0 || dk >= oz)
            voro_fatal_error("Constructing periodic image for nonexistent point", VOROPP_INTERNAL_ERROR);
        if (dk >= ez && dk < wz)
        {
            if (dj < ey || dj >= wy)
                create_side_image(di, dj, dk);
        }
        else
            create_vertical_image(di, dj, dk);
    }
    void create_side_image(int di, int dj, int dk);
    void create_vertical_image(int di, int dj, int dk);
    void put_image(int reg, int fijk, int l, double dx, double dy, double dz);
    inline void remap(int& ai, int& aj, int& ak, int& ci, int& cj, int& ck, double& x, double& y, double& z, int& ijk);
};

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_C_LOOP_BASE_HH
#define VOROPP_C_LOOP_BASE_HH

#include "config.hh"
namespace voro
{

/** \brief Base class for looping over particles in a container.
 *
 * This class forms the base of all classes that can loop over a subset of
 * particles in a contaner in some order. When initialized, it stores constants
 * about the corresponding container geometry. It also contains a number of
 * routines for interrogating which particle currently being considered by the
 * loop, which are common between all of the derived classes. */
class c_loop_base
{
public:
    /** The number of blocks in the x direction. */
    const int nx;
    /** The number of blocks in the y direction. */
    const int ny;
    /** The number of blocks in the z direction. */
    const int nz;
    /** A constant, set to the value of nx multiplied by ny, which
     * is used in the routines that step through blocks in
     * sequence. */
    const int nxy;
    /** A constant, set to the value of nx*ny*nz, which is used in
     * the routines that step through blocks in sequence. */
    const int nxyz;
    /** The number of floating point numbers per particle in the
     * associated container data structure. */
    const int ps;
    /** A pointer to the particle position information in the
     * associated container data structure. */
    double** p;
    /** A pointer to the particle ID information in the associated
     * container data structure. */
    int** id;
    /** A pointer to the particle counts in the associated
     * container data structure. */
    int* co;
    /** The current x-index of the block under consideration by the
     * loop. */
    int i;
    /** The current y-index of the block under consideration by the
     * loop. */
    int j;
    /** The current z-index of the block under consideration by the
     * loop. */
    int k;
    /** The current index of the block under consideration by the
     * loop. */
    int ijk;
    /** The index of the particle under consideration within the current
     * block. */
    int q;
    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    template<class c_class>
    c_loop_base(c_class& con)
        : nx(con.nx)
        , ny(con.ny)
        , nz(con.nz)
        , nxy(con.nxy)
        , nxyz(con.nxyz)
        , ps(con.ps)
        , p(con.p)
        , id(con.id)
        , co(con.co)
    {
    }
    /** Returns the position vector of the particle currently being
     * considered by the loop.
     * \param[out] (x,y,z) the position vector of the particle. */
    inline void pos(double& x, double& y, double& z)
    {
        double* pp = p[ijk] + ps * q;
        x = *(pp++);
        y = *(pp++);
        z = *pp;
    }
    /** Returns the ID, position vector, and radius of the particle
     * currently being considered by the loop.
     * \param[out] pid the particle ID.
     * \param[out] (x,y,z) the position vector of the particle.
     * \param[out] r the radius of the particle. If no radius
     * 		 information is available the default radius
     * 		 value is returned. */
    inline void pos(int& pid, double& x, double& y, double& z, double& r)
    {
        pid = id[ijk][q];
        double* pp = p[ijk] + ps * q;
        x = *(pp++);
        y = *(pp++);
        z = *pp;
        r = ps == 3 ? default_radius : *(++pp);
    }
    /** Returns the x position of the particle currently being
     * considered by the loop. */
    inline double x()
    {
        return p[ijk][ps * q];
    }
    /** Returns the y position of the particle currently being
     * considered by the loop. */
    inline double y()
    {
        return p[ijk][ps * q + 1];
    }
    /** Returns the z position of the particle currently being
     * considered by the loop. */
    inline double z()
    {
        return p[ijk][ps * q + 2];
    }
    /** Returns the ID of the particle currently being considered
     * by the loop. */
    inline int pid()
    {
        return id[ijk][q];
    }
};

} // namespace voro

#endif

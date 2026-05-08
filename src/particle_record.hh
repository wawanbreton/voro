// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_PARTICLE_RECORD_HH
#define VOROPP_PARTICLE_RECORD_HH

namespace voro
{

/** \brief Structure for holding information about a particle.
 *
 * This small structure holds information about a single particle, and is used
 * by several of the routines in the voro_compute template for passing
 * information by reference between functions. */
struct particle_record
{
    /** The index of the block that the particle is within. */
    int ijk;
    /** The number of particle within its block. */
    int l;
    /** The x-index of the block. */
    int di;
    /** The y-index of the block. */
    int dj;
    /** The z-index of the block. */
    int dk;
};

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_C_LOOP_SUBSET_MODE_HH
#define VOROPP_C_LOOP_SUBSET_MODE_HH

namespace voro
{

/** A type associated with a c_loop_subset class, determining what type of
 * geometrical region to loop over. */
enum c_loop_subset_mode
{
    sphere,
    box,
    no_check
};

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_WALL_HH
#define VOROPP_WALL_HH

#include "voronoicell.hh"
#include "voronoicell_neighbor.hh"
namespace voro
{

/** \brief Pure virtual class from which wall objects are derived.
 *
 * This is a pure virtual class for a generic wall object. A wall object
 * can be specified by deriving a new class from this and specifying the
 * functions.*/
class wall
{
public:
    virtual ~wall()
    {
    }
    /** A pure virtual function for testing whether a point is
     * inside the wall object. */
    virtual bool point_inside(double x, double y, double z) = 0;
    /** A pure virtual function for cutting a cell without
     * neighbor-tracking with a wall. */
    virtual bool cut_cell(voronoicell& c, double x, double y, double z) = 0;
    /** A pure virtual function for cutting a cell with
     * neighbor-tracking enabled with a wall. */
    virtual bool cut_cell(voronoicell_neighbor& c, double x, double y, double z) = 0;
};

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_WALL_LIST_HH
#define VOROPP_WALL_LIST_HH

#include "wall.hh"
namespace voro
{

/** \brief A class for storing a list of pointers to walls.
 *
 * This class stores a list of pointers to wall classes. It contains several
 * simple routines that make use of the wall classes (such as telling whether a
 * given position is inside all of the walls or not). It can be used by itself,
 * but also forms part of container_base, for associating walls with this
 * class. */
class wall_list
{
public:
    /** An array holding pointers to wall objects. */
    wall** walls;
    /** A pointer to the next free position to add a wall pointer.
     */
    wall** wep;
    wall_list();
    ~wall_list();
    /** Adds a wall to the list.
     * \param[in] w the wall to add. */
    inline void add_wall(wall* w)
    {
        if (wep == wel)
            increase_wall_memory();
        *(wep++) = w;
    }
    /** Adds a wall to the list.
     * \param[in] w a reference to the wall to add. */
    inline void add_wall(wall& w)
    {
        add_wall(&w);
    }
    void add_wall(wall_list& wl);
    /** Determines whether a given position is inside all of the
     * walls on the list.
     * \param[in] (x,y,z) the position to test.
     * \return True if it is inside, false if it is outside. */
    inline bool point_inside_walls(double x, double y, double z)
    {
        for (wall** wp = walls; wp < wep; wp++)
            if (! ((*wp)->point_inside(x, y, z)))
                return false;
        return true;
    }
    /** Cuts a Voronoi cell by all of the walls currently on
     * the list.
     * \param[in] c a reference to the Voronoi cell class.
     * \param[in] (x,y,z) the position of the cell.
     * \return True if the cell still exists, false if the cell is
     * deleted. */
    template<class c_class>
    bool apply_walls(c_class& c, double x, double y, double z)
    {
        for (wall** wp = walls; wp < wep; wp++)
            if (! ((*wp)->cut_cell(c, x, y, z)))
                return false;
        return true;
    }
    void deallocate();

protected:
    void increase_wall_memory();
    /** A pointer to the limit of the walls array, used to
     * determine when array is full. */
    wall** wel;
    /** The current amount of memory allocated for walls. */
    int current_wall_size;
};

} // namespace voro

#endif

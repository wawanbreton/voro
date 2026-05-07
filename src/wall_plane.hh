// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_WALL_PLANE_HH
#define VOROPP_WALL_PLANE_HH

#include "wall.hh"
namespace voro
{

/** \brief A class representing a plane wall object.
 *
 * This class represents a single plane wall object. */
struct wall_plane : public wall
{
public:
    /** Constructs a plane wall object.
     * \param[in] (xc_,yc_,zc_) a normal vector to the plane.
     * \param[in] ac_ a displacement along the normal vector.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    wall_plane(double xc_, double yc_, double zc_, double ac_, int w_id_ = -99)
        : w_id(w_id_)
        , xc(xc_)
        , yc(yc_)
        , zc(zc_)
        , ac(ac_)
    {
    }
    bool point_inside(double x, double y, double z);
    template<class v_cell>
    bool cut_cell_base(v_cell& c, double x, double y, double z);
    bool cut_cell(voronoicell& c, double x, double y, double z)
    {
        return cut_cell_base(c, x, y, z);
    }
    bool cut_cell(voronoicell_neighbor& c, double x, double y, double z)
    {
        return cut_cell_base(c, x, y, z);
    }

private:
    const int w_id;
    const double xc, yc, zc, ac;
};

} // namespace voro

#endif

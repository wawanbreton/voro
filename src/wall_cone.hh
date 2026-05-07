// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_WALL_CONE_HH
#define VOROPP_WALL_CONE_HH

#include <cmath>
#include "wall.hh"
namespace voro
{

/** \brief A class representing a conical wall object.
 *
 * This class represents a cone wall object. */
struct wall_cone : public wall
{
public:
    /** Constructs a cone wall object.
     * \param[in] (xc_,yc_,zc_) the apex of the cone.
     * \param[in] (xa_,ya_,za_) a vector pointing along the axis of
     *			    the cone.
     * \param[in] ang the angle (in radians) of the cone, measured
     *		  from the axis.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    wall_cone(double xc_, double yc_, double zc_, double xa_, double ya_, double za_, double ang, int w_id_ = -99)
        : w_id(w_id_)
        , xc(xc_)
        , yc(yc_)
        , zc(zc_)
        , xa(xa_)
        , ya(ya_)
        , za(za_)
        , asi(1 / (xa_ * xa_ + ya_ * ya_ + za_ * za_))
        , gra(tan(ang))
        , sang(sin(ang))
        , cang(cos(ang))
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
    const double xc, yc, zc, xa, ya, za, asi, gra, sang, cang;
};

} // namespace voro

#endif

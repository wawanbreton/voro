// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_WALL_CYLINDER_HH
#define VOROPP_WALL_CYLINDER_HH

#include "wall.hh"
namespace voro
{

/** \brief A class representing a cylindrical wall object.
 *
 * This class represents a open cylinder wall object. */
struct wall_cylinder : public wall
{
public:
    /** Constructs a cylinder wall object.
     * \param[in] (xc_,yc_,zc_) a point on the axis of the
     *			    cylinder.
     * \param[in] (xa_,ya_,za_) a vector pointing along the
     *			    direction of the cylinder.
     * \param[in] rc_ the radius of the cylinder
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    wall_cylinder(double xc_, double yc_, double zc_, double xa_, double ya_, double za_, double rc_, int w_id_ = -99)
        : w_id(w_id_)
        , xc(xc_)
        , yc(yc_)
        , zc(zc_)
        , xa(xa_)
        , ya(ya_)
        , za(za_)
        , asi(1 / (xa_ * xa_ + ya_ * ya_ + za_ * za_))
        , rc(rc_)
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
    const double xc, yc, zc, xa, ya, za, asi, rc;
};

} // namespace voro

#endif

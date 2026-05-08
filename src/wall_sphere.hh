// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_WALL_SPHERE_HH
#define VOROPP_WALL_SPHERE_HH

#include "wall.hh"
namespace voro
{

/** \brief A class representing a spherical wall object.
 *
 * This class represents a spherical wall object. */
struct wall_sphere : public wall
{
public:
    /** Constructs a spherical wall object.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking.
     * \param[in] (xc_,yc_,zc_) a position vector for the sphere's
     * 			    center.
     * \param[in] rc_ the radius of the sphere. */
    wall_sphere(double xc_, double yc_, double zc_, double rc_, int w_id_ = -99)
        : w_id(w_id_)
        , xc(xc_)
        , yc(yc_)
        , zc(zc_)
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
    const double xc, yc, zc, rc;
};

} // namespace voro

#endif

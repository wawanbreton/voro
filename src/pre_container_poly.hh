// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_PRE_CONTAINER_POLY_HH
#define VOROPP_PRE_CONTAINER_POLY_HH

#include <cstdio>
#include "common.hh"
#include "particle_order.hh"
#include "container_poly.hh"
#include "pre_container_base.hh"
namespace voro
{

/** \brief A class for storing an arbitrary number of particles with radius
 * information, prior to setting up a container geometry.
 *
 * The pre_container_poly class is an extension of the pre_container_base class
 * for cases when particle radius information is available. */
class pre_container_poly : public pre_container_base
{
public:
    /** The class constructor sets up the geometry of container,
     * initializing the minimum and maximum coordinates in each
     * direction.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
     *                                                container is periodic in
     *                                                each coordinate direction. */
    pre_container_poly(
        double ax_,
        double bx_,
        double ay_,
        double by_,
        double az_,
        double bz_,
        bool xperiodic_,
        bool yperiodic_,
        bool zperiodic_)
        : pre_container_base(ax_, bx_, ay_, by_, az_, bz_, xperiodic_, yperiodic_, zperiodic_, 4) {};
    void put(int n, double x, double y, double z, double r);
    void import(FILE* fp = stdin);
    /** Imports particles from a file.
     * \param[in] filename the name of the file to read from. */
    inline void import(const char* filename)
    {
        FILE* fp = safe_fopen(filename, "r");
        import(fp);
        fclose(fp);
    }
    void setup(container_poly& con);
    void setup(particle_order& vo, container_poly& con);
};

} // namespace voro

#endif

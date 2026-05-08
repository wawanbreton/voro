// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_C_LOOP_SUBSET_HH
#define VOROPP_C_LOOP_SUBSET_HH

#include "c_loop_base.hh"
#include "c_loop_subset_mode.hh"
namespace voro
{

/** \brief Class for looping over a subset of particles in a container.
 *
 * This class can loop over a subset of particles in a certain geometrical
 * region within the container. The class can be set up to loop over a
 * rectangular box or sphere. It can also rectangular group of internal
 * computational blocks. */
class c_loop_subset : public c_loop_base
{
public:
    /** The current mode of operation, determining whether tests
     * should be applied to particles to ensure they are within a
     * certain geometrical object. */
    c_loop_subset_mode mode;
    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    template<class c_class>
    c_loop_subset(c_class& con)
        : c_loop_base(con)
        , ax(con.ax)
        , ay(con.ay)
        , az(con.az)
        , sx(con.bx - ax)
        , sy(con.by - ay)
        , sz(con.bz - az)
        , xsp(con.xsp)
        , ysp(con.ysp)
        , zsp(con.zsp)
        , xperiodic(con.xperiodic)
        , yperiodic(con.yperiodic)
        , zperiodic(con.zperiodic)
    {
    }
    void setup_sphere(double vx, double vy, double vz, double r, bool bounds_test = true);
    void setup_box(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, bool bounds_test = true);
    void setup_intbox(int ai_, int bi_, int aj_, int bj_, int ak_, int bk_);
    bool start();
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc()
    {
        do
        {
            q++;
            while (q >= co[ijk])
            {
                q = 0;
                if (! next_block())
                    return false;
            }
        } while (mode != no_check && out_of_bounds());
        return true;
    }

private:
    const double ax, ay, az, sx, sy, sz, xsp, ysp, zsp;
    const bool xperiodic, yperiodic, zperiodic;
    double px, py, pz, apx, apy, apz;
    double v0, v1, v2, v3, v4, v5;
    int ai, bi, aj, bj, ak, bk;
    int ci, cj, ck, di, dj, dk, inc1, inc2;
    inline int step_mod(int a, int b)
    {
        return a >= 0 ? a % b : b - 1 - (b - 1 - a) % b;
    }
    inline int step_div(int a, int b)
    {
        return a >= 0 ? a / b : -1 + (a + 1) / b;
    }
    inline int step_int(double a)
    {
        return a < 0 ? int(a) - 1 : int(a);
    }
    void setup_common();
    bool next_block();
    bool out_of_bounds();
};

} // namespace voro

#endif

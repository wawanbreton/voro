// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_C_LOOP_ALL_PERIODIC_HH
#define VOROPP_C_LOOP_ALL_PERIODIC_HH

#include "c_loop_base.hh"
namespace voro
{

/** \brief A class for looping over all particles in a container_periodic or
 * container_periodic_poly class.
 *
 * Since the container_periodic and container_periodic_poly classes have a
 * fundamentally different memory organization, the regular loop classes cannot
 * be used with them. */
class c_loop_all_periodic : public c_loop_base
{
public:
    /** The constructor copies several necessary constants from the
     * base periodic container class.
     * \param[in] con the periodic container class to use. */
    template<class c_class>
    c_loop_all_periodic(c_class& con)
        : c_loop_base(con)
        , ey(con.ey)
        , ez(con.ez)
        , wy(con.wy)
        , wz(con.wz)
        , ijk0(nx * (ey + con.oy * ez))
        , inc2(2 * nx * con.ey + 1)
    {
    }
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    inline bool start()
    {
        i = 0;
        j = ey;
        k = ez;
        ijk = ijk0;
        q = 0;
        while (co[ijk] == 0)
            if (! next_block())
                return false;
        return true;
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc()
    {
        q++;
        if (q >= co[ijk])
        {
            q = 0;
            do
            {
                if (! next_block())
                    return false;
            } while (co[ijk] == 0);
        }
        return true;
    }

private:
    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    int ey;
    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    int ez;
    /** The upper y index (exclusive) of the primary domain within
     * the block structure. */
    int wy;
    /** The upper z index (exclusive) of the primary domain within
     * the block structure. */
    int wz;
    /** The index of the (0,0,0) block within the block structure.
     */
    int ijk0;
    /** A value to increase ijk by when the z index is increased.
     */
    int inc2;
    /** Updates the internal variables to find the next
     * computational block with any particles.
     * \return True if another block is found, false if there are
     * no more blocks. */
    inline bool next_block()
    {
        i++;
        if (i == nx)
        {
            i = 0;
            j++;
            if (j == wy)
            {
                j = ey;
                k++;
                if (k == wz)
                    return false;
                ijk += inc2;
            }
            else
                ijk++;
        }
        else
            ijk++;
        return true;
    }
};

} // namespace voro

#endif

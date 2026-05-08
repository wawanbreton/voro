// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#ifndef VOROPP_C_LOOP_ORDER_HH
#define VOROPP_C_LOOP_ORDER_HH

#include "c_loop_base.hh"
#include "particle_order.hh"
namespace voro
{

/** \brief Class for looping over all of the particles specified in a
 * pre-assembled particle_order class.
 *
 * The particle_order class can be used to create a specific order of particles
 * within the container. This class can then loop over these particles in this
 * order. The class is particularly useful in cases where the ordering of the
 * output must match the ordering of particles as they were inserted into the
 * container. */
class c_loop_order : public c_loop_base
{
public:
    /** A reference to the ordering class to use. */
    particle_order& vo;
    /** A pointer to the current position in the ordering class. */
    int* cp;
    /** A pointer to the end position in the ordering class. */
    int* op;
    /** The constructor copies several necessary constants from the
     * base class, and sets up a reference to the ordering class to
     * use.
     * \param[in] con the container class to use.
     * \param[in] vo_ the ordering class to use. */
    template<class c_class>
    c_loop_order(c_class& con, particle_order& vo_)
        : c_loop_base(con)
        , vo(vo_)
        , nx(con.nx)
        , nxy(con.nxy)
    {
    }
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    inline bool start()
    {
        cp = vo.o;
        op = vo.op;
        if (cp != op)
        {
            ijk = *(cp++);
            decode();
            q = *(cp++);
            return true;
        }
        else
            return false;
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc()
    {
        if (cp == op)
            return false;
        ijk = *(cp++);
        decode();
        q = *(cp++);
        return true;
    }

private:
    /** The number of computational blocks in the x direction. */
    const int nx;
    /** The number of computational blocks in a z-slice. */
    const int nxy;
    /** Takes the current block index and computes indices in the
     * x, y, and z directions. */
    inline void decode()
    {
        k = ijk / nxy;
        int ijkt = ijk - nxy * k;
        j = ijkt / nx;
        i = ijkt - j * nx;
    }
};

} // namespace voro

#endif

#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <projamrex.H>

using namespace amrex;

void Projamrex::initBCs()
{
    /* ====== BC TYPES ======
        0: Periodic
        1: Slip walls
        2: No-slip walls
    */

    auto setBCType = [this] (int bcint, Orientation ori)
    {
        m_bc_den[ori] = 1.0;
        m_bc_press[ori] = 1.0;
        m_bc_vel[ori][0] = 0.0;
        m_bc_vel[ori][1] = 0.0;
        m_bc_trac[ori].resize(m_ntrac, 0.0);

        if(bcint == 1)
        {
            m_bc_type[ori] = BC::slipWall;
        }
        else if(bcint == 2)
        {
            m_bc_type[ori] = BC::noSlipWall;

        }
        else
        {
            m_bc_type[ori] = BC::undefined;
        }

        if(geom.isPeriodic(ori.coordDir()))
        {
            if(m_bc_type[ori] == BC::undefined)
            {
                m_bc_type[ori] == BC::periodic;
            }
            else
            {
                amrex::Abort("wrong BC type for periodic geom");
            }
        }
    };
    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    ParmParse ppBC("init");
    ppBC.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    ppBC.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);

    setBCType(lo_bc[0], Orientation(Direction::x, Orientation::low));
    setBCType(lo_bc[1], Orientation(Direction::y, Orientation::low));
    setBCType(hi_bc[0], Orientation(Direction::x, Orientation::high));
    setBCType(hi_bc[1], Orientation(Direction::y, Orientation::high));


    /* =================== VELOCITY =================== */
    m_bcrec_vel.resize(AMREX_SPACEDIM);
    for(OrientationIter oit; oit; oit++)
    {
        Orientation ori = oit();
        int dir = ori.coordDir();
        Orientation::Side side = ori.faceDir();
        auto const bct = m_bc_type[ori];
        if(bct == BC::slipWall)
        {
            if(side == Orientation::low)
            {
                // Tangential directions have hoextrap
                m_bcrec_vel[0].setLo(dir, BCType::hoextrap);
                m_bcrec_vel[1].setLo(dir, BCType::hoextrap);

                // Only normal direction has ext_dir
                m_bcrec_vel[dir].setLo(dir, BCType::ext_dir);
            }
            else
            {
                   // Tangential directions have hoextrap
                   m_bcrec_vel[0].setHi(dir, BCType::hoextrap);
                   m_bcrec_vel[1].setHi(dir, BCType::hoextrap);

                   // Only normal direction has ext_dir
                   m_bcrec_vel[dir].setHi(dir, BCType::ext_dir);
              }
        }
        else if(bct == BC::noSlipWall)
        {
            if(side == Orientation::low)
            {
                m_bcrec_vel[0].setLo(dir, BCType::ext_dir);
                m_bcrec_vel[1].setLo(dir, BCType::ext_dir);
            }
            else
            {
                m_bcrec_vel[0].setHi(dir, BCType::ext_dir);
                m_bcrec_vel[1].setHi(dir, BCType::ext_dir);
            }
        }
        else if(bct == BC::periodic)
        {
            if(side == Orientation::low)
            {
                m_bcrec_vel[0].setLo(dir, BCType::int_dir);
                m_bcrec_vel[1].setLo(dir, BCType::int_dir);
            }
            else
            {
                m_bcrec_vel[0].setHi(dir, BCType::int_dir);
                m_bcrec_vel[1].setHi(dir, BCType::int_dir);
            }
        }
    }


    /* =================== DENSITY =================== */

    m_bcrec_den.resize(1);
    for(OrientationIter oit; oit; oit++)
    {
        Orientation ori = oit();
        int dir = ori.coordDir();
        Orientation::Side side = ori.faceDir();
        auto const bct = m_bc_type[ori];
        if(bct == BC::slipWall)
        {
            if(side == Orientation::low)
            {
                m_bcrec_den[0].setLo(dir, BCType::hoextrap);

            }
            else
            {
                m_bcrec_den[0].setHi(dir, BCType::hoextrap);
            }
        }
        else if(bct == BC::noSlipWall)
        {
            if(side == Orientation::low)
            {
                m_bcrec_den[0].setLo(dir, BCType::foextrap);
            }
            else
            {
                m_bcrec_den[0].setLo(dir, BCType::foextrap);
            }
        }
        else if(bct == BC::periodic)
        {
            if(side == Orientation::low)
            {
                m_bcrec_den[0].setLo(dir, BCType::int_dir);
            }
            else
            {
                m_bcrec_den[0].setLo(dir, BCType::int_dir);
            }
        }
    }

    /* =================== PRESSURE =================== */

    m_bcrec_press.resize(1);
    for(OrientationIter oit; oit; oit++)
    {
        Orientation ori = oit();
        int dir = ori.coordDir();
        Orientation::Side side = ori.faceDir();
        auto const bct = m_bc_type[ori];
        if(bct == BC::slipWall)
        {
            if(side == Orientation::low)
            {
                m_bcrec_press[0].setLo(dir, BCType::foextrap);

            }
            else
            {
                m_bcrec_press[0].setHi(dir, BCType::foextrap);
            }
        }
        else if(bct == BC::noSlipWall)
        {
            if(side == Orientation::low)
            {
                m_bcrec_press[0].setLo(dir, BCType::foextrap);
            }
            else
            {
                m_bcrec_press[0].setLo(dir, BCType::foextrap);
            }
        }
        else if(bct == BC::periodic)
        {
            if(side == Orientation::low)
            {
                m_bcrec_press[0].setLo(dir, BCType::int_dir);
            }
            else
            {
                m_bcrec_press[0].setLo(dir, BCType::int_dir);
            }
        }
    }


    // =================== TRACER ===================


    if(m_ntrac > 0)
    {
        m_bcrec_trac.resize(m_ntrac);
        for (OrientationIter oit; oit; ++oit)
        {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = m_bc_type[ori];
            if(bct == BC::slipWall)
            {
                if(side == Orientation::low)
                {
                    for (auto& b : m_bcrec_trac) b.setLo(dir, BCType::hoextrap);
                }
                else
                {
                    for (auto& b : m_bcrec_trac) b.setHi(dir, BCType::hoextrap);
                }
            }
            else if(bct == BC::noSlipWall)
            {
                if(side == Orientation::low)
                {
                    for (auto& b : m_bcrec_trac) b.setLo(dir, BCType::foextrap);
                }
                else
                {
                    for (auto& b : m_bcrec_trac) b.setHi(dir, BCType::foextrap);
                }
            }
            else if(bct == BC::periodic)
            {
                if(side == Orientation::low)
                {
                    for (auto& b : m_bcrec_trac) b.setLo(dir, BCType::int_dir);
                }
                else
                {
                    for (auto& b : m_bcrec_trac) b.setHi(dir, BCType::int_dir);
                }
            }
        }
    }




}

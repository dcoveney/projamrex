#include <projamrex.H>
#include <bcFuncs.H>

using namespace amrex;

void Projamrex::fillphysbc_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    PhysBCFunct<GpuBndryFuncFab<ProjVelFill> > physbc(geom, m_bcrec_vel,
                                                        ProjVelFill(testNumber, m_bc_vel));
    physbc.FillBoundary(vel, 0, AMREX_SPACEDIM, IntVect(ng), time, 0);
}

#include <AMReX_BC_TYPES.H>
#include <projamrex.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MacProjector.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_AmrCore.H>

using namespace amrex;


Array<amrex::LinOpBCType,AMREX_SPACEDIM>
Projamrex::fetchBCsForMacProjection (Orientation::Side side) const noexcept
{
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
    {
    if (geom.isPeriodic(dim))
    {
        r[dim] = LinOpBCType::Periodic;
    }
    }
    return r;
}


void
Projamrex::doMacProjection(Array<MultiFab, AMREX_SPACEDIM>& vel)
{
    int finest_level = parent->finestLevel();

    Array<MultiFab, AMREX_SPACEDIM> beta;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 3);
        beta[idim].setVal(1.0);  // set beta to 1
    }

    Vector<Array<MultiFab*,AMREX_SPACEDIM> > mac_vec(finest_level+1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        mac_vec[lev][0] = &vel[0];
        mac_vec[lev][1] = &vel[1];
        mac_vec[lev][2] = &vel[2];
    }


    LPInfo lp_info;

    MacProjector macproj(mac_vec,       // face-based velocity
                         {amrex::GetArrOfConstPtrs(beta)}, // beta
                         {geom},                           // the geometry object
                         lp_info);

     auto lo_bc = fetchBCsForMacProjection(Orientation::low);
     auto hi_bc = fetchBCsForMacProjection(Orientation::high);

     macproj.setDomainBC(lo_bc, hi_bc);
     macproj.setVerbose(2);

     // Define the relative tolerance
     Real reltol = 1.e-8;

     // Define the absolute tolerance; note that this argument is optional
     Real abstol = 1.e-15;

     macproj.project(reltol,abstol);


}

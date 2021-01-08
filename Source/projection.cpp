#include <AMReX_BC_TYPES.H>
#include <projamrex.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_NodalProjector.H>
#include <AMReX_MLABecLaplacian.H>
using namespace amrex;

Array<amrex::LinOpBCType,AMREX_SPACEDIM>
Projamrex::fetchBCsForProjection (Orientation::Side side) const noexcept
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
Projamrex::doNodalProjection(MultiFab& vel, MultiFab& pressure, MultiFab& S_cc, MultiFab& rho, Real dt)
{
    // Solve elliptic equation for the nodal pressure field
    // D(grad(p)/rho) = D(U*)/dt = S_cc/dt
    // D(U*) is nodal, vel is cell-centred
    // The projection is performed via vel = vel - grad(p)/rho

    LPInfo lp_info;

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (parent->Geom(0).isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            if (phys_bc.lo(idim) == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc.lo(idim) == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            if (phys_bc.hi(idim) == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc.hi(idim) == Inflow) {
                mlmg_hibc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
        }
    }

    Real dt_inv = 1.0/dt;

    MLNodeLaplacian matrix({geom}, {grids}, {dmap}, lp_info);
    matrix.setDomainBC(mlmg_lobc, mlmg_hibc);
    lp_info.setMaxCoarseningLevel(30);
    matrix.setGaussSeidel(true);
    matrix.setHarmonicAverage(false);
    S_cc.mult(dt_inv,1);
    const BoxArray& ndgrids = amrex::convert(grids, IntVect{1,1});
    MultiFab S_nd(ndgrids, dmap, 1, 1);
    S_nd.setVal(0.0); // Set it to zero for this example
    MultiFab rhs(ndgrids, dmap, 1, 1);
    matrix.compRHS({&rhs}, {&vel}, {&S_nd}, {&S_cc});
    MultiFab sigma(grids, dmap, 1, 1);
    sigma.setVal(1.0);

    // matrix.setSigma(0, sigma);
    matrix.setSigma(0, rho);

    MLMG nodal_solver(matrix);
    // nodal_solver.setMaxIter(mg_maxiter);
    // nodal_solver.setBottomMaxIter(mg_bottom_maxiter);

    nodal_solver.setVerbose(2);
    Real reltol = 1.0e-10;
    Real abstol = 1.0e-15;
    nodal_solver.solve({&pressure}, {&rhs}, reltol, abstol);
    MultiFab projFluxes(grids, dmap, AMREX_SPACEDIM, 1);
    projFluxes.setVal(0.0);
    nodal_solver.getFluxes({&projFluxes});
    if(projFluxes.contains_nan()) Abort("-1/rho * gp contains nan");

    projFluxes.mult(dt, 1);
    MultiFab::Add(vel, projFluxes, 0, 0, AMREX_SPACEDIM, 0);
}

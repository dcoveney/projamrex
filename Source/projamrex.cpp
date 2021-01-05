
#include <projamrex.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <PROB_NS_F.H>
#include <NS_BC.H>
#include <AMReX_Math.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_NodalProjector.H>

using namespace amrex;

int      Projamrex::verbose         = 1;
Real     Projamrex::cfl             = 0.9;
int      Projamrex::do_reflux       = 0;

int      Projamrex::NUM_STATE       = 2;  // Four variables in the state
int      Projamrex::NUM_GROW        = 4;  // number of ghost cells
int      Projamrex::num_state_type  = 2;  // number of state types: state, pressure
Real     Projamrex::mu              = 0.0;
int      Projamrex::m_ntrac         = 1;
int      Projamrex::testNumber      = 0;
BCRec    Projamrex::phys_bc;


std::string probType                = "";
int    Nstate                          = 5;
int    Npress                          = 1;
int    Density                         = 0;
int    Xvel                            = 1;
int    Yvel                            = 2;
int    Trac                            = 3;
int    magvort                        = 4;
int    Press                           = 0;



#ifdef AMREX_PARTICLES
std::unique_ptr<AmrTracerParticleContainer> Projamrex::TracerPC =  nullptr;
int Projamrex::do_tracers                       =  0;
#endif

//
//Default constructor.  Builds invalid object.
//
Projamrex::Projamrex ()
{
    flux_reg = 0;

}

//
//The basic constructor.
//
Projamrex::Projamrex (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{

    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);

}

//
//The destructor.
//
Projamrex::~Projamrex ()
{
    delete flux_reg;
}


static
void
set_x_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM == 3)
static
void
set_z_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}
#endif

static
void
set_scalar_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_temp_bc (BCRec&       bc,
             const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,temp_bc[lo_bc[i]]);
        bc.setHi(i,temp_bc[hi_bc[i]]);
    }
}

static
void
set_pressure_bc (BCRec&       bc,
                 const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,press_bc[lo_bc[i]]);
        bc.setHi(i,press_bc[hi_bc[i]]);
    }
}

static
void
set_divu_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,divu_bc[lo_bc[i]]);
        bc.setHi(i,divu_bc[hi_bc[i]]);
    }
}

static
void
set_dsdt_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,dsdt_bc[lo_bc[i]]);
        bc.setHi(i,dsdt_bc[hi_bc[i]]);
    }
}

//
//Restart from a checkpoint file.
//
void
Projamrex::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

void
Projamrex::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
#ifdef AMREX_PARTICLES
  if (do_tracers and level == 0) {
    TracerPC->Checkpoint(dir, "Tracer", true);
  }
#endif
}

//
//Write a plotfile to specified directory.
//
void
Projamrex::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                            VisMF::How         how)
{

    AmrLevel::writePlotFile (dir,os,how);

#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      TracerPC->Checkpoint(dir, "Tracer", true);
    }
#endif
}

//
//Define data descriptors.
//

typedef StateDescriptor::BndryFunc BndryFunc;

void
Projamrex::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    // read_params();
    ParmParse pp_init("init");
    pp_init.query("testNumber", testNumber);
    pp_init.query("vel_visc_coef", mu);
    pp_init.query("probType", probType);
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,Nstate,
			   &cell_cons_interp);
    desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                           StateDescriptor::Interval,1,1,
			   &node_bilinear_interp,true);

   Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
   pp_init.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
   pp_init.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);
   for (int i = 0; i < AMREX_SPACEDIM; i++)
   {
       phys_bc.setLo(i,lo_bc[i]);
       phys_bc.setHi(i,hi_bc[i]);
   }

        BCRec bc;

        set_scalar_bc(bc,phys_bc);


        desc_lst.setComponent(State_Type, Density, "density", bc,
                  BndryFunc(FORT_DENFILL));

        set_x_vel_bc(bc,phys_bc);

        desc_lst.setComponent(State_Type, Xvel, "vel_x", bc,
                  BndryFunc(FORT_XVELFILL));

        set_y_vel_bc(bc,phys_bc);
        desc_lst.setComponent(State_Type, Yvel, "vel_y", bc,
                  BndryFunc(FORT_YVELFILL));

         set_scalar_bc(bc,phys_bc);

        desc_lst.setComponent(State_Type, Trac, "tracer", bc,
                  BndryFunc(FORT_ADVFILL));
      desc_lst.setComponent(State_Type, magvort, "mag_vort", bc,
                BndryFunc(FORT_ADVFILL));
        set_pressure_bc(bc,phys_bc);

        desc_lst.setComponent(Press_Type, Press, "pressure",bc,
                  BndryFunc(FORT_PRESFILL));



}

//
//Cleanup data descriptors at end of run.
//
void
Projamrex::variableCleanUp ()
{
    desc_lst.clear();
#ifdef AMREX_PARTICLES
    TracerPC.reset();
#endif
}

//
//Initialize grid data at problem start-up.
//
void
Projamrex::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //

    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);
    Real cur_time   = state[State_Type].curTime();

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        RealBox    gridloc = RealBox(box,geom.CellSize(),geom.ProbLo());
        FArrayBox& Sfab = S_new[mfi];
        initialState(Sfab, box, gridloc, dx, testNumber);

          // initdata(&level, &cur_time, AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
		  //  BL_TO_FORTRAN_3D(S_new[mfi]), &NUM_STATE, AMREX_ZFILL(dx),
		  //  AMREX_ZFILL(prob_lo));
    }

    // don't need pressure until the projection
    P_new.setVal(0.0);

#ifdef AMREX_PARTICLES
    init_particles();
#endif
IntVect cell(AMREX_D_DECL(16,32,0));
print_state(S_new, cell);
    if (verbose) {
	amrex::Print() << "Done initializing the level " << level
                       << " data " << std::endl;
    }
}

//
//Initialize data on this level from another Projamrex (during regrid).
//
void
Projamrex::init (AmrLevel &old)
{
    Projamrex* oldlev = (Projamrex*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);

    FillPatch(old, S_new, 0, cur_time, State_Type, 0, S_new.nGrow());
}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
Projamrex::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    // MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, S_new.nGrow());
}

//
//Advance grids at this level in time.
//
Real
Projamrex::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{
    MultiFab& S_mm = get_new_data(State_Type);

    std::ofstream amrexEulerOut;
    for (int k = 0; k < num_state_type; k++)
    {
	bool has_old_data = state[k].hasOldData();
        state[k].allocOldData();
        // state[k].swapTimeLevels(dt);
    }

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_old = get_old_data(Press_Type);
    MultiFab& P_new = get_new_data(Press_Type);

    IntVect cell(AMREX_D_DECL(16,32,0));
    Print() << "New data:\n";
    print_state(S_new, cell);
    Print() << "Old data:\n";
    print_state(S_old, cell);
    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time = state[State_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;

    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
	fine = &getFluxReg(level+1);
	fine->setVal(0.0);
    }

    if (do_reflux && level > 0) {
	current = &getFluxReg(level);
    }

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux)
    {
	for (int j = 0; j < BL_SPACEDIM; j++)
	{
	    BoxArray ba = S_new.boxArray();
	    ba.surroundingNodes(j);
	    fluxes[j].define(ba, dmap, NUM_STATE, 0);
	}
    }

    int ngmac = 2;
    int ngrow = 3;
    // State with ghost cells

    MultiFab velCC(grids, dmap, AMREX_SPACEDIM, 1);


    MultiFab Sborder(grids, dmap, Nstate, ngrow);
    // MultiFab::Copy(Sborder, S_new, 0, 0, Nstate, S_new.nGrow());
    FillPatch(*this, Sborder, ngrow, cur_time, State_Type, 0, Nstate);

    MultiFab rho_inv(grids, dmap, 1, 1);
    MultiFab divu_cc(grids, dmap, 1, 1);

    // MF to hold the mac velocity
    MultiFab Umac[BL_SPACEDIM], uface[BL_SPACEDIM], flux[BL_SPACEDIM];
    Array<MultiFab, AMREX_SPACEDIM> uAdvN, uAdvT;
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
      const BoxArray& edgeba = getEdgeBoxArray(i);
      Umac[i].define(convert(grids,IntVect::TheDimensionVector(i)), dmap, 1, ngmac);
      uAdvN[i].define(convert(grids,IntVect::TheDimensionVector(i)), dmap, 1, ngmac);
      uface[i].define(convert(grids,IntVect::TheDimensionVector(i)), dmap, 1, ngmac);
      flux[i].define(convert(grids,IntVect::TheDimensionVector(i)), dmap, Nstate, ngmac);
      Umac[i].setVal(0.0);
      uAdvN[i].setVal(0.0);
      uface[i].setVal(0.0);
    }

    const BoxArray& edgebay = getEdgeBoxArray(1);
    const BoxArray& edgebax = getEdgeBoxArray(0);
    uAdvT[0].define(convert(grids,IntVect::TheDimensionVector(1)), dmap, 1, ngmac);
    uAdvT[1].define(convert(grids,IntVect::TheDimensionVector(0)), dmap, 1, ngmac);
    uAdvT[0].setVal(0.0);
    uAdvT[1].setVal(0.0);
    //Defining face-centred MultiFabs for the staggered velocity field
    // average_cellcenter_to_face(GetArrOfPtrs(faceVel[1]),velCC_y, geom);
    // average_cellcenter_to_face(GetArrOfPtrs(faceVel[1]),velCC[1], geom);

    if(S_new.contains_nan())
    {
        Abort("Before advection: S_new contains nan!");
    }
    if(Sborder.contains_nan())
    {
        Abort("Before advection: Sborder contains nan!");
    }
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
        const Box& bx = mfi.tilebox();

	    const FArrayBox& statein = Sborder[mfi];
	    FArrayBox& stateout      =   S_new[mfi];
        FArrayBox& rhoinvFAB     = rho_inv[mfi];

        setFaceVelocity(uface[0][mfi], uface[1][mfi], statein, bx);


        computeNormalAdvectionVelocity(statein, uface[0][mfi], uface[1][mfi],
            uAdvN[0][mfi], uAdvN[1][mfi], bx, dx, dt);

    }
    }

    MultiFab::Copy(Sborder, S_new, 0, 0, Nstate, S_new.nGrow());
    FillPatch(*this, Sborder, ngrow, cur_time, State_Type, 0, Nstate);

    doMacProjection(uAdvN);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
        const Box& bx = mfi.tilebox();

        const FArrayBox& statein = Sborder[mfi];
        FArrayBox& stateout      =   S_new[mfi];

        advectScalars(statein, stateout, uface[0][mfi], uface[1][mfi], uAdvN[0][mfi], uAdvN[1][mfi],
            flux[0][mfi], flux[1][mfi], bx, dx, dt);
        // scalarAdvection(statein, stateout, flux[0][mfi], flux[1][mfi], uface[0][mfi], uface[1][mfi], bx, dx, dt);
        computeTangentialAdvectionVelocity(statein, uAdvN[0][mfi], uAdvN[1][mfi],
            uAdvT[0][mfi], uAdvT[1][mfi], bx, dx, dt);
        // velocityAdvection(statein, stateout, flux[0][mfi], flux[1][mfi],uface[0][mfi], uface[1][mfi], bx, dx, dt);
        if (do_reflux) {
		for (int i = 0; i < BL_SPACEDIM ; i++)
		    fluxes[i][mfi].copy<RunOn::Host>(flux[i][mfi],mfi.nodaltilebox(i));
	    }
    }
    }
    if(uAdvN[0].contains_nan() || uAdvN[1].contains_nan())
    {
        Abort("advection velocity contains nan");
    }
    if(uAdvT[0].contains_nan() || uAdvT[1].contains_nan())
    {
        Abort("advection velocity contains nan");
    }

    MultiFab::Copy(Sborder, S_new, 0, 0, Nstate, S_new.nGrow());
    FillPatch(*this, Sborder, ngrow, cur_time, State_Type, 0, Nstate);

    MultiFab uAdvCC(grids, dmap, AMREX_SPACEDIM, 0);
    average_face_to_cellcenter(uAdvCC, 0,
        Array<MultiFab const*,AMREX_SPACEDIM>{{AMREX_D_DECL(&uAdvN[0],&uAdvN[1],&uAdvN[2])}});

    if(uAdvCC.contains_nan())
    {
        Abort("uAdvCC contains nan!");
    }
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
        {

    	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    	{
            const FArrayBox& statein = Sborder[mfi];
            FArrayBox& stateout      =   S_new[mfi];
            const Box& bx = mfi.tilebox();

            advectVelocity(uAdvCC[mfi], uAdvN[0][mfi], uAdvN[1][mfi],
                uAdvT[0][mfi], uAdvT[1][mfi], statein, stateout, bx, dx, dt);
        }
        }


        MultiFab::Copy(Sborder, S_new, 0, 0, Nstate, S_new.nGrow());
        FillPatch(*this, Sborder, ngrow, cur_time, State_Type, 0, Nstate);
        if(S_new.contains_nan())
        {
            Abort("S_new contains nan!");
        }
        if(Sborder.contains_nan())
        {
            Abort("Sborder contains nan!");
        }
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
        {
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    	{
            const FArrayBox& statein = Sborder[mfi];
            FArrayBox& stateout      =   S_new[mfi];
            FArrayBox& rhoinvFAB     = rho_inv[mfi];
            const Box& bx = mfi.tilebox();

             velocityDiffusion(statein, stateout, bx, dx, dt);

            Dim3 lo = lbound(bx);
            Dim3 hi = ubound(bx);
            Array4<Real const> const& stateArray = statein.array();
            Array4<Real> const& rhoinvArr  = rhoinvFAB.array();
            for(int k = lo.z; k <= hi.z; k++)
            {
                for(int j = lo.y-1; j <= hi.y+1; j++)
                {
                    for(int i = lo.x-1; i <= hi.x+1; i++)
                    {
                        rhoinvArr(i,j,k,0) = 1.0/(stateArray(i,j,k,Density));

                    }
                }
            }
        }
        }

        if(S_new.contains_nan())
        {
            Abort("After diffusion: S_new contains nan!");
        }
        if(Sborder.contains_nan())
        {
            Abort("After diffusion: Sborder contains nan!");
        }
    // found u_tmp, projection stuff goes here
    // need to compute a node-centred approximation to the divergence of u_tmp
    // average u_tmp to CCs, use central differencing for gradient, average to nodes
    // then apply the projection

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
        {

    	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    	{
            const Box& bx = mfi.tilebox();

            const FArrayBox& statein = Sborder[mfi];
            FArrayBox& divuFAB       = divu_cc[mfi];
            computeDivergence(divuFAB, statein, bx, dx);

        }
        }
    MultiFab::Copy(Sborder, S_new, 0, 0, Nstate, S_new.nGrow());
    FillPatch(*this, Sborder, ngrow, cur_time, State_Type, 0, Nstate);

    MultiFab::Copy(velCC, Sborder, Xvel,0,AMREX_SPACEDIM,1);
    if(S_new.contains_nan())
    {
        Abort("Before projection: S_new contains nan!");
    }
    if(Sborder.contains_nan())
    {
        Abort("Before projection: Sborder contains nan!");
    }

    MultiFab rho(grids, dmap, 1, S_new.nGrow());
    MultiFab::Copy(rho, S_new, Density,0,1,S_new.nGrow());
    Print() << "BEFORE PROJECTION\n";
    Print() << "DIVU\n";
    print_state(divu_cc, cell);
    Print() << "VEL\n";
    print_state(velCC, cell);
    if(velCC.contains_nan())
    {
        Abort("Before projection: velCC contains nan!");
    }
    doNodalProjection(velCC,P_new,divu_cc,rho_inv,dt);
    //
    if(velCC.contains_nan())
    {
        Abort("After projection: velCC contains nan!");
    }
    MultiFab::Copy(S_new, velCC, 0, Xvel, AMREX_SPACEDIM, S_new.nGrow());
    MultiFab::Copy(Sborder, S_new, 0, 0, Nstate, S_new.nGrow());
    FillPatch(*this, Sborder, ngrow, cur_time, State_Type, 0, Nstate);

    if(S_new.contains_nan())
    {
        Abort("After projection: S_new contains nan!");
    }
    if(Sborder.contains_nan())
    {
        Abort("After projection: Sborder contains nan!");
    }
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
        {
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    	{
            const FArrayBox& statein = Sborder[mfi];
            FArrayBox& derout      =   S_new[mfi];
            const Box& bx = mfi.tilebox();
            computeMagVorticity(statein, derout, bx, dx);


        }
        }
    if (do_reflux) {
    if (current) {
        for (int i = 0; i < BL_SPACEDIM ; i++)
        current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    }
    if (fine) {
        for (int i = 0; i < BL_SPACEDIM ; i++)
        fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    }
    }
    // MultiFab SBorder(grids, dmap, NUM_STATE, NUM_GROW);



// #ifdef AMREX_PARTICLES
//     if (TracerPC) {
//       TracerPC->AdvectWithUmac(Umac, level, dt);
//     }
// #endif
    // amrex::Abort();
    return dt;

}

void
Projamrex::initialState(FArrayBox& statein, const Box& bx,
    RealBox gridloc, const Real* dx, int testNumber)
{
    Array4<Real> const& stateArray = statein.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    const Real* xlo = gridloc.lo();
    const Real* xhi = gridloc.hi();
    ParmParse pp("init");
    Real x, y, z;
    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {
                x = xlo[0] + dx[0]*((i-lo.x) + 0.5);
                y = xlo[1] + dx[1]*((j-lo.y) + 0.5);
                #if(AMREX_SPACEDIM==3)
                z = xlo[2] + dx[2]*((k-lo.z) + 0.5);
                #endif

                if(testNumber == 1)
                {
                    stateArray(i,j,k,Xvel) = 1.0;
                    stateArray(i,j,k,Yvel) = 0.0;
                    // stateArray(i,j,k,Density) = 1.0 - 0.5*tanh((x-0.5)/0.01);
                    if(x < 0.5)
                    {
                        stateArray(i,j,k,Density) = 1.0;
                        stateArray(i,j,k,Trac) = 1.0;
                    }
                    else
                    {
                        stateArray(i,j,k,Density) = 0.5;
                        stateArray(i,j,k,Trac) = 0.0;
                    }
                }
                else if(testNumber == 2)
                {
                    stateArray(i,j,k, Xvel) = tanh(30.0*(0.25-fabs(y-0.5)));
                    stateArray(i,j,k, Yvel) = 0.05*sin(2*M_PI*x);
                    stateArray(i,j,k, Density) = 1.0;
                    stateArray(i,j,k,Trac) = 0.0;
                    stateArray(i,j,k,magvort) = 0.0;

                }
                else if(testNumber == 3)
                {
                    stateArray(i,j,k,Xvel) = 0.0;
                    stateArray(i,j,k,Yvel) = 0.0;
                    stateArray(i,j,k,Density) = 1.0;
                    stateArray(i,j,k, Trac) = 0.0;
                    stateArray(i,j,k, magvort) = 0.0;
                }
            }
        }
    }
}

void
Projamrex::setFaceVelocity(FArrayBox& vel_x, FArrayBox& vel_y, const FArrayBox& statein, const Box& bx)
{
    Array4<Real> const& u = vel_x.array();
    Array4<Real> const& v = vel_y.array();
    Array4<Real const> const& stateArr = statein.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y-1; j <= hi.y+1; j++)
        {
            for(int i = lo.x-2; i <= hi.x+1; i++)
            {
                if(isnan(stateArr(i,j,k,Xvel)))
                {
                    Abort("Xvel is nan");
                }
                u(i+1,j,k) = 0.5*(stateArr(i+1,j,k,Xvel)+stateArr(i,j,k,Xvel));

            }
        }
    }
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y-2; j <= hi.y+1; j++)
        {
            for(int i = lo.x-1; i <= hi.x+1; i++)
            {
                v(i,j+1,k) = 0.5*(stateArr(i,j+1,k,Yvel)+stateArr(i,j,k,Yvel));
            }
        }
    }
}

void Projamrex::advectScalars(const FArrayBox& statein, FArrayBox& stateout, const FArrayBox& u_edge,
const FArrayBox& v_edge, const FArrayBox& u_half, const FArrayBox& v_half, FArrayBox& flux_x, FArrayBox& flux_y,
    const Box& bx, const Real* dx, const Real dt)
 {
     // CC density
     Array4<Real const> const& den = statein.array();
     Array4<Real> const& den_out = stateout.array();
     // x component of normal advective velocity on x-faces
     Array4<Real const> const& u_hf = u_half.array();
     // y component of normal advective velocity on y-faces
     Array4<Real const> const& v_hf = v_half.array();
    // x velocity on x faces
     Array4<Real const> const& u_ed = u_edge.array();
     // y velocity on y faces
     Array4<Real const> const& v_ed = v_edge.array();
     Array4<Real> const& flx = flux_x.array();
     Array4<Real> const& fly = flux_y.array();

     Dim3 lo = lbound(bx);
     Dim3 hi = ubound(bx);

     Real sLnph_x, sRnph_x, u_adv, unorm, vnorm, sLnph_y, sRnph_y, snph_x, snph_y;
     Vector<Real> sx(2), sty(2), vtan(2);
     Vector<Real> sy(2), stx(2), utan(2);

     for(int k = lo.z; k <= hi.z; k++)
     {
         for(int j = lo.y-1; j <= hi.y+1; j++)
         {
             for(int i = lo.x-2; i <= hi.x+1; i++)
             {
                 unorm = 0.5*(u_hf(i+1,j,k)+u_hf(i,j,k));
                 sx[0] = 0.5*(den(i+1,j,k,Density)-den(i-1,j,k,Density));
                 vtan[0] = 0.5*(v_hf(i,j+1,k)+v_hf(i,j,k));
                 sty[0] = (vtan[0] > 0.0) ? (den(i,j,k,Density)-den(i,j-1,k,Density)) :
                             (den(i,j+1,k,Density)-den(i,j,k,Density));
                 sLnph_x = den(i,j,k,Density) + fmin(0.5*(1.0-unorm*(dt/dx[0])),0.5)*sx[0]
                         - (dt/(2.0*dx[0]))*sty[0]*vtan[0];

                 sx[1] = 0.5*(den(i+2,j,k,Density)-den(i,j,k,Density));
                 vtan[1] = 0.5*(v_hf(i+1,j+1,k)+v_hf(i+1,j,k));
                 sty[1] = (vtan[1] > 0.0) ? (den(i+1,j,k,Density)-den(i+1,j-1,k,Density)) :
                             (den(i+1,j+1,k,Density)-den(i+1,j,k,Density));
                 sRnph_x = den(i,j,k,Density) + fmax(0.5*(-1.0-unorm*(dt/dx[0])),-0.5)*sx[1]
                         - (dt/(2.0*dx[0]))*sty[1]*vtan[1];

                snph_x = upwind(sLnph_x, sRnph_x, u_ed(i+1,j,k));
                flx(i+1,j,k,Density) = u_hf(i+1,j,k)*snph_x;
            }
         }
     }
     for(int k = lo.z; k <= hi.z; k++)
     {
         for(int j = lo.y-2; j <= hi.y+1; j++)
         {
             for(int i = lo.x-1; i <= hi.x+1; i++)
             {
                 vnorm = 0.5*(v_hf(i,j+1,k)+v_hf(i,j,k));
                 sy[0] = 0.5*(den(i,j+1,k,Density)-den(i,j-1,k,Density));
                 utan[0] = 0.5*(u_hf(i+1,j,k)+u_hf(i,j,k));
                 stx[0] = (utan[0] > 0.0) ? (den(i,j,k,Density)-den(i-1,j,k,Density)) :
                             (den(i+1,j,k,Density)-den(i,j,k,Density));
                 sLnph_y = den(i,j,k,Density) + fmin(0.5*(1.0-vnorm*(dt/dx[1])),0.5)*sy[0]
                         - (dt/(2.0*dx[1]))*stx[0]*utan[0];

                 sy[1] = 0.5*(den(i,j+2,k,Density)-den(i,j,k,Density));
                 utan[1] = 0.5*(u_hf(i+1,j+1,k)+u_hf(i,j+1,k));
                 stx[1] = (utan[1] > 0.0) ? (den(i,j+1,k,Density)-den(i-1,j+1,k,Density)) :
                             (den(i+1,j+1,k,Density)-den(i,j+1,k,Density));
                 sRnph_y = den(i,j,k,Density) + fmax(0.5*(-1.0-vnorm*(dt/dx[1])),-0.5)*sy[1]
                         - (dt/(2.0*dx[1]))*stx[1]*utan[1];

                snph_y = upwind(sLnph_y, sRnph_y, v_ed(i,j+1,k));
                fly(i,j+1,k,Density) = v_hf(i,j+1,k)*snph_y;
            }
         }
     }
     for(int k = lo.z; k <= hi.z; k++)
     {
         for(int j = lo.y; j <= hi.y; j++)
         {
             for(int i = lo.x; i <= hi.x; i++)
             {
                 if(den_out(i,j,k,Density)==0.0)
                 {
                     Abort("before scalar advection update: zero density");
                 }

                 den_out(i,j,k,Density) = den(i,j,k,Density) ;//+
                       // (dt/dx[0])*(flx(i,j,k, Density)-flx(i+1,j,k, Density)) +
                       // (dt/dx[1])*(fly(i,j,k, Density)-fly(i,j+1,k, Density));

                 if(den_out(i,j,k,Density)==0.0)
                 {
                     Abort("after scalar advection update: zero density");
                 }

             }
         }
     }

 }

void
Projamrex::scalarAdvection(const FArrayBox& statein, FArrayBox& stateout, FArrayBox& flux_x, FArrayBox& flux_y,
 const FArrayBox& vel_x, const FArrayBox& vel_y, const Box& bx, const Real* dx, const Real dt)
{
    Array4<Real const> const& den = statein.array();
    Array4<Real const> const& u = vel_x.array();
    Array4<Real const> const& v = vel_y.array();
    Array4<Real> const& flx = flux_x.array();
    Array4<Real> const& fly = flux_y.array();
    Array4<Real> const& stateArr = stateout.array();

    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x+1; i++)
            {
                if(u(i,j,k) > 0)
                {
                    flx(i,j,k,0) = u(i,j,k)*den(i-1,j,k,Density);
                }
                else
                {
                    flx(i,j,k,0) = u(i,j,k)*den(i,j,k,Density);
                }
            }
        }
    }
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y+1; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                if(v(i,j,k) > 0)
                {
                    fly(i,j,k,0) = v(i,j,k)*den(i,j-1,k,Density);
                }
                else
                {
                    fly(i,j,k,0) = v(i,j,k)*den(i,j,k,Density);
                }
            }
        }
    }

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                if(stateArr(i,j,k,Density)==0.0)
                {
                    Abort("before scalar advection update: zero density");
                }

                stateArr(i,j,k,Density) = den(i,j,k,Density) +
                    (dt/dx[0])*(flx(i,j,k,0)-flx(i+1,j,k,0)) +
                    (dt/dx[1])*(fly(i,j,k,0)-fly(i,j+1,k,0));

                if(stateArr(i,j,k,Density)==0.0)
                {
                    Abort("after scalar advection update: zero density");
                }

            }
        }
    }
}

void
Projamrex::computeNormalAdvectionVelocity(const FArrayBox& statein, const FArrayBox& vel_x, const FArrayBox& vel_y,
        FArrayBox& advVelx, FArrayBox& advVely, const Box& bx, const Real* dx, const Real dt)
{
    // CC vel
    Array4<Real const> const& vel = statein.array();
    // u on x-faces
    Array4<Real const> const& u_mac = vel_x.array();
    // v on y-faces
    Array4<Real const> const& v_mac = vel_y.array();
    // face-centred normal advection velocities, to store u_{i+1/2,j}^{n+1/2}, v_{i,j+1,2}^{n+1/2}
    Array4<Real> const& adv_x = advVelx.array();
    Array4<Real> const& adv_y = advVely.array();

    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    // index 0 = LEFT STATES, 1 = RIGHT STATES
    Vector<Real> ux(2), uty(2), vtan(2);
    Vector<Real> vy(2), vtx(2), utan(2);
    Real uLnph, uRnph, u_adv, unorm, vnorm, vLnph, vRnph;

    // x velocities
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y-1; j <= hi.y+1; j++)
        {
            for(int i = lo.x-2; i <= hi.x+1; i++)
            {
                unorm = 0.5*(u_mac(i+1,j,k)+u_mac(i,j,k));
                ux[0] = 0.5*(vel(i+1,j,k,Xvel)-vel(i-1,j,k,Xvel));
                vtan[0] = 0.5*(v_mac(i,j+1,k)+v_mac(i,j,k));
                uty[0] = (vtan[0] > 0.0) ? (vel(i,j,k,Xvel)-vel(i,j-1,k,Xvel)) :
                            (vel(i,j+1,k,Xvel)-vel(i,j,k,Xvel));
                if(isnan(ux[0]))
                {
                    Abort("uxL is nan");
                }
                if(isnan(unorm))
                {
                    Abort("unorn is nan");
                }
                uLnph = vel(i,j,k,Xvel) + fmin(0.5*(1.0-unorm*(dt/dx[0])),0.5)*ux[0]
                        - (dt/(2.0*dx[0]))*uty[0]*vtan[0];

                if(isnan(uLnph))
                {
                    Abort("left extrap state is nan");
                }
                ux[1] = 0.5*(vel(i+2,j,k,Xvel)-vel(i,j,k,Xvel));
                vtan[1] = 0.5*(v_mac(i+1,j+1,k)+v_mac(i+1,j,k));
                uty[1] = (vtan[1] > 0.0) ? (vel(i+1,j,k,Xvel)-vel(i+1,j-1,k,Xvel)) :
                            (vel(i+1,j+1,k,Xvel)-vel(i+1,j,k,Xvel));

                if(isnan(ux[1]))
                {
                    Abort("uxR is nan");
                }
                if(isnan(uRnph))
                {
                    Abort("right extrap state is nan");
                }
                uRnph = vel(i+1,j,k,Xvel) + fmax(0.5*(-1.0-unorm*(dt/dx[0])),-0.5)*ux[1]
                        - (dt/(2.0*dx[0]))*uty[1]*vtan[1];

                adv_x(i+1,j,k) = upwind(uLnph, uRnph, u_mac(i+1,j,k));

            }
        }
    }

    //y velocities

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y-2; j <= hi.y+1; j++)
        {
            for(int i = lo.x-1; i <= hi.x+1; i++)
            {
                // left edge
                vnorm = 0.5*(v_mac(i,j+1,k)+v_mac(i,j,k));
                vy[0] = 0.5*(vel(i,j+1,k,Yvel)-vel(i,j-1,k,Yvel));
                utan[0] = 0.5*(u_mac(i+1,j,k)+u_mac(i,j,k));
                vtx[0] = (utan[0] > 0.0) ? (vel(i,j,k,Yvel)-vel(i-1,j,k,Yvel)) :
                            (vel(i+1,j,k,Yvel)-vel(i,j,k,Yvel));

                vLnph = vel(i,j,k,Yvel) + fmin(0.5*(1.0-vnorm*(dt/dx[1])),0.5)*vy[0]
                        - (dt/(2.0*dx[1]))*vtx[0]*utan[0];

                if(isnan(vy[0]))
                {
                    Abort("vyL is nan");
                }
                if(isnan(vnorm))
                {
                    Abort("vnorn is nan");
                }
                if(isnan(vLnph))
                {
                    Abort("v left extrap state is nan");
                }
                // right edge
                vy[1] = 0.5*(vel(i,j+2,k,Yvel)-vel(i,j,k,Yvel));
                utan[1] = 0.5*(u_mac(i+1,j+1,k)+u_mac(i,j+1,k));
                vtx[1] = (utan[1] > 0.0) ? (vel(i,j+1,k,Yvel)-vel(i-1,j+1,k,Yvel)) :
                            (vel(i+1,j+1,k,Yvel)-vel(i,j+1,k,Yvel));

                vRnph = vel(i,j+1,k,Yvel) + fmax(0.5*(-1.0-vnorm*(dt/dx[1])),-0.5)*vy[1]
                        - (dt/(2.0*dx[1]))*vtx[1]*utan[1];

                if(isnan(vy[1]))
                {
                    Abort("vyR is nan");
                }
                if(isnan(vRnph))
                {
                    Abort("v right extrap state is nan");
                }

                adv_y(i,j+1,k) = upwind(vLnph, vRnph, v_mac(i,j+1,k));

            }
        }
    }


}

void
Projamrex::computeTangentialAdvectionVelocity(const FArrayBox& statein, const FArrayBox& vel_x, const FArrayBox& vel_y,
        FArrayBox& advVelx, FArrayBox& advVely, const Box& bx, const Real* dx, const Real dt)
{
    // CC vel
    Array4<Real const> const& vel = statein.array();
    // x component of normal advective velocity on x-faces
    Array4<Real const> const& u_half = vel_x.array();
    // y component of normal advective velocity on y-faces
    Array4<Real const> const& v_half = vel_y.array();
    // face-centred tangential advection velocities, to store u_{i,j+1/2}^{n+1/2}, v_{i+1/2,j}^{n+1/2}
    Array4<Real> const& adv_x = advVelx.array();
    Array4<Real> const& adv_y = advVely.array();

    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    // index 0 = LEFT STATES, 1 = RIGHT STATES
    Vector<Real> uy(2), utx(2), utan(2);
    Vector<Real> vx(2), vty(2), vtan(2);
    Real uLnph, uRnph, u_adv, unorm, vnorm, vLnph, vRnph;


    // compute v_{i+1/2,j}^{n+1/2}
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y-1; j <= hi.y+1; j++)
        {
            for(int i = lo.x-2; i <= hi.x+1; i++)
            {
                unorm = 0.5*(u_half(i+1,j,k)+u_half(i,j,k));
                vx[0] = 0.5*(vel(i+1,j,k,Yvel)-vel(i-1,j,k,Yvel));
                vtan[0] = 0.5*(v_half(i,j+1,k)+v_half(i,j,k));
                vty[0] = (vtan[0] > 0.0) ? (vel(i,j,k,Yvel)-vel(i,j-1,k,Yvel)) :
                            (vel(i,j+1,k,Yvel)-vel(i,j,k,Yvel));

                vLnph = vel(i,j,k,Yvel) + fmin(0.5*(1.0-unorm*(dt/dx[0])),0.5)*vx[0]
                        - (dt/(2.0*dx[0]))*vty[0]*vtan[0];

                vx[1] = 0.5*(vel(i+2,j,k,Yvel)-vel(i,j,k,Yvel));
                vtan[1] = 0.5*(v_half(i+1,j+1,k)+v_half(i+1,j,k));
                vty[1] = (vtan[1] > 0.0) ? (vel(i+1,j,k,Yvel)-vel(i+1,j-1,k,Yvel)) :
                            (vel(i+1,j+1,k,Yvel)-vel(i+1,j,k,Yvel));

                vRnph = vel(i+1,j,k,Yvel) + fmax(0.5*(-1.0-unorm*(dt/dx[0])),-0.5)*vx[1]
                        - (dt/(2.0*dx[0]))*vty[1]*vtan[1];

                adv_y(i+1,j,k) = upwind(vLnph, vRnph, u_half(i+1,j,k));
            }
        }
    }

    // compute u_{i,j+1/2}^{n+1/2}

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y-2; j <= hi.y+1; j++)
        {
            for(int i = lo.x-1; i <= hi.x+1; i++)
            {
                vnorm = 0.5*(v_half(i,j+1,k)+v_half(i,j,k));
                uy[0] = 0.5*(vel(i,j+1,k,Xvel)-vel(i,j-1,k,Xvel));
                utan[0] = 0.5*(u_half(i+1,j,k)+u_half(i,j,k));
                utx[0] = (utan[0] > 0.0) ? (vel(i,j,k,Xvel)-vel(i-1,j,k,Xvel)) :
                            (vel(i+1,j,k,Xvel)-vel(i,j,k,Xvel));

                uLnph = vel(i,j,k,Xvel) + fmin(0.5*(1.0-vnorm*(dt/dx[1])),0.5)*uy[0]
                        - (dt/(2.0*dx[1]))*utx[0]*utan[0];

                vnorm = 0.5*(v_half(i,j+1,k)+v_half(i,j,k));
                uy[1] = 0.5*(vel(i,j+2,k,Xvel)-vel(i,j,k,Xvel));
                utan[1] = 0.5*(u_half(i+1,j+1,k)+u_half(i,j+1,k));
                utx[1] = (utan[1] > 0.0) ? (vel(i,j+1,k,Xvel)-vel(i-1,j+1,k,Xvel)) :
                            (vel(i+1,j+1,k,Xvel)-vel(i,j+1,k,Xvel));

                uRnph = vel(i,j+1,k,Xvel) + fmax(0.5*(-1.0-vnorm*(dt/dx[1])),-0.5)*uy[1]
                        - (dt/(2.0*dx[1]))*utx[1]*utan[1];

                adv_x(i,j+1,k) = upwind(uLnph, uRnph, v_half(i,j+1,k));

            }
        }
    }


}

void
Projamrex::advectVelocity(const FArrayBox& advVelCC, const FArrayBox& u_halfx, const FArrayBox& v_halfy,
    const FArrayBox& u_halfy, const FArrayBox& v_halfx, const FArrayBox& statein, FArrayBox& stateout,
     const Box& bx, const Real* dx, const Real dt)
{
    // CC vel
    Array4<Real const> const& advVel = advVelCC.array();
    // u_{i+1/2,j}^half
    Array4<Real const> const& u_halfN = u_halfx.array();
    // v_{i,j+1/2}^half
    Array4<Real const> const& v_halfN = v_halfy.array();
    // u_{i,j+1/2}^half
    Array4<Real const> const& u_halfT = u_halfy.array();
    // v_{i+1/2,j}^half
    Array4<Real const> const& v_halfT = v_halfx.array();
    // to store u*
    Array4<Real> const& vel_out           = stateout.array();
    Array4<Real const> const& vel_in  = statein.array();
    //to store [(u dot grad)u]_{i,j}^{n+1/2} [(u dot grad)v]_{i,j}^{n+1/2}
    double advx, advy;

    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                advx = advVel(i,j,k,0)*(u_halfN(i+1,j,k)-u_halfN(i,j,k))/(dx[0]) +
                    advVel(i,j,k,1)*(u_halfT(i,j+1,k)-u_halfT(i,j,k))/(dx[1]);
                advy = advVel(i,j,k,0)*(v_halfT(i+1,j,k)-v_halfT(i,j,k))/(dx[0]) +
                    advVel(i,j,k,1)*(v_halfN(i,j+1,k)-v_halfN(i,j,k))/(dx[1]);

                vel_out(i,j,k,Xvel) = vel_in(i,j,k,Xvel) - dt*advx;
                vel_out(i,j,k,Yvel) = vel_in(i,j,k,Yvel) - dt*advy;

            }
        }
    }
}

void
Projamrex::velocityAdvection(const FArrayBox& statein, FArrayBox& stateout, FArrayBox& flux_x, FArrayBox& flux_y,
 const FArrayBox& vel_x, const FArrayBox& vel_y, const Box& bx, const Real* dx, const Real dt)
{
    Array4<Real const> const& vel = statein.array();
    Array4<Real const> const& u_mac = vel_x.array();
    Array4<Real const> const& v_mac = vel_y.array();
    Array4<Real> const& flx = flux_x.array();
    Array4<Real> const& fly = flux_y.array();
    Array4<Real> const& stateArr = stateout.array();

    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    // index 0 = LEFT STATES, 1 = RIGHT STATES
    Vector<Real> unorm(2), ux(2), uyt(2), vtan(2);
    Real uLnph, uRnph, u_adv;
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x-1; i <= hi.x; i++)
            {
                if(u_mac(i,j,k) > 0)
                {
                    flx(i,j,k,Xvel) = u_mac(i,j,k)*vel(i-1,j,k,Xvel);
                    flx(i,j,k,Yvel) = u_mac(i,j,k)*vel(i,j,k,Yvel);
                }
                else
                {
                    flx(i,j,k,Xvel) = u_mac(i,j,k)*vel(i,j,k,Xvel);
                    flx(i,j,k,Yvel) = u_mac(i,j,k)*vel(i,j,k,Yvel);
                }

            }
        }
    }
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y+1; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {

                if(v_mac(i,j,k) > 0)
                {
                    fly(i,j,k,Yvel) = v_mac(i,j,k)*vel(i,j-1,k,Yvel);
                    fly(i,j,k,Xvel) = v_mac(i,j,k)*vel(i,j,k,Xvel);
                }
                else
                {
                    fly(i,j,k,Yvel) = v_mac(i,j,k)*vel(i,j,k,Yvel);
                    fly(i,j,k,Xvel) = v_mac(i,j,k)*vel(i,j,k,Xvel);
                }
            }
        }
    }

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                for(int n = Xvel; n <= Yvel; n++)
                {
                    stateArr(i,j,k,n) = vel(i,j,k,n) +
                         (dt/dx[0])*(flx(i,j,k,n)-flx(i+1,j,k,n)) +
                         (dt/dx[1])*(fly(i,j,k,n)-fly(i,j+1,k,n));

                }

            }
        }
    }
}



// void Projamrex::advectionENO(const FArrayBox& statein, FArrayBox& stateout,const Box& bx, const Real* dx, const Real dt)
// {
//
// }

void
Projamrex::velocityDiffusion(const FArrayBox& statein, FArrayBox& stateout,
    const Box& bx, const Real* dx, const Real dt)
{
    Array4<Real const> const& vel = statein.array();
    Array4<Real> const& stateArr = stateout.array();
    double rho, uxx, uyy, vxx, vyy;
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                uxx = (vel(i+1,j,k,Xvel)-2.0*vel(i,j,k,Xvel)+vel(i-1,j,k,Xvel))/(dx[0]*dx[0]);
                uyy = (vel(i,j+1,k,Xvel)-2.0*vel(i,j,k,Xvel)+vel(i,j-1,k,Xvel))/(dx[1]*dx[1]);
                vxx = (vel(i+1,j,k,Yvel)-2.0*vel(i,j,k,Yvel)+vel(i-1,j,k,Yvel))/(dx[0]*dx[0]);
                vyy = (vel(i,j+1,k,Yvel)-2.0*vel(i,j,k,Yvel)+vel(i,j-1,k,Yvel))/(dx[1]*dx[1]);
                rho = vel(i,j,k,Density);

                // if(rho==0.0)
                // {
                //     Abort("zero density");
                // }
                stateArr(i,j,k,Xvel) = vel(i,j,k,Xvel) + (dt/rho)*mu*(uxx+uyy);
                stateArr(i,j,k,Yvel) = vel(i,j,k,Yvel) + (dt/rho)*mu*(vxx+vyy);


            }
        }
    }
}

void
Projamrex::computeAdvectionTerm(FArrayBox& vel_x, FArrayBox& vel_y, FArrayBox& adv_x, FArrayBox& adv_y,
        const Box& bx, const Real* dx)
{
    Array4<Real const> const& velArr_x = vel_x.array();
    Array4<Real const> const& velArr_y = vel_y.array();
    Array4<Real> const& advArr_x = adv_x.array();
    Array4<Real> const& advArr_y = adv_y.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    Real uC, uL, uR, uU, uD, vR, vC, vRD, uLU, vL, vU, vD;

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                uC = velArr_x(i,j,k,0); uR = velArr_x(i+1,j,k,0);
                uL = velArr_x(i-1,j,k,0); uU = velArr_x(i,j+1,k,0);
                uD = velArr_x(i,j-1,k,0); vC = velArr_y(i, j, k, 0);
                vR = velArr_y(i+1,j, k, 0); vL = velArr_y(i-1,j,k,0);
                vU = velArr_y(i,j+1, k, 0); vD = velArr_y(i,j-1,k,0);
                uLU = velArr_x(i-1,j+1, k, 0); vRD = velArr_y(i+1,j-1,k,0);

                advArr_x(i,j,k,0) = (1.0/dx[0])*(0.25*(uR+uC)*(uR+uC)-0.25*(uC+uL)*(uC+uL))
                    + (1.0/dx[1])*(0.25*(uU+uC)*(vR+vC)-0.25*(uC+uD)*(vRD+vC));
                advArr_y(i,j,k,0) = (1.0/dx[0])*(0.25*(uC+uU)*(vC+vR)-0.25*(uLU+uL)*(vC+vL))
                    + (1.0/dx[1])*(0.25*(vR+vC)*(vR+vC)-0.25*(vC+vD)*(vC+vD));
            }
        }
    }


}
void
Projamrex::computeTemporaryVelocity(FArrayBox& adv_x, FArrayBox& adv_y, FArrayBox& tempvel_x, FArrayBox& tempvel_y,
        FArrayBox& vel_x, FArrayBox& vel_y,const Box& bx, const Real* dx, const Real dt)
{
    Array4<Real> const& tempvelArr_x = tempvel_x.array();
    Array4<Real> const& tempvelArr_y = tempvel_y.array();
    Array4<Real const> const& velArr_x = vel_x.array();
    Array4<Real const> const& velArr_y = vel_y.array();
    Array4<Real const> const& advArr_x = adv_x.array();
    Array4<Real const> const& advArr_y = adv_y.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    Real uC, uL, uR, uU, uD, vR, vC, vRD, uLU, vL, vU, vD;

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                tempvelArr_x(i,j,k,0) = velArr_x(i,j,k,0) + dt*(-advArr_x(i,j,k,0));
                tempvelArr_y(i,j,k,0) = velArr_y(i,j,k,0) + dt*(-advArr_y(i,j,k,0));
            }
        }
    }

}

void
Projamrex::computeAdvectiveFlux(const FArrayBox& vel, FArrayBox& rho,
             FArrayBox& flux, const Box& bx, const Real* dx, int dim)
{
    Array4<Real const> const& velArr = vel.array();
    Array4<Real const> const& den = rho.array();
    Array4<Real> const& fluxArr = flux.array();

    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    double v = 0.0;
    if(dim == 0)
    {

        for(int k = lo.z; k <= hi.z; ++k)
        {
            for(int j = lo.y; j <= hi.y; ++j)
            {
                for(int i = lo.x; i <= hi.x+1; ++i)
                {
                    v = velArr(i,j,k,0);
                    if(v > 0)
                    {
                        fluxArr(i,j,k,0) = v*den(i-1,j,k,Density);
                    }
                    else
                    {
                        fluxArr(i,j,k,0) = v*den(i,j,k,Density);
                    }
                }
            }
        }
    }
    else if(dim == 1)
    {
        for(int k = lo.z; k <= hi.z; ++k)
        {
            for(int j = lo.y; j <= hi.y+1; ++j)
            {
                for(int i = lo.x; i <= hi.x; ++i)
                {
                    v = velArr(i,j,k,1);
                    if(v > 0)
                    {
                        fluxArr(i,j,k,0) = 1.0*den(i,j-1,k, Density);
                    }
                    else
                    {
                        fluxArr(i,j,k,0) = v*den(i,j,k, Density);
                    }
                }
            }
        }
    }
}

void
Projamrex::computeDivergence(FArrayBox& divu, const FArrayBox& vel, const Box& bx, const Real* dx)
{
    Array4<Real const> const& velArr = vel.array();
    Array4<Real> const& divuArr = divu.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    double dudx, dvdy;
    for(int k = lo.z; k <= hi.z; ++k)
    {
        for(int j = lo.y-1; j <= hi.y+1; ++j)
        {
            for(int i = lo.x-1; i <= hi.x+1; ++i)
            {
                dudx = (velArr(i+1, j, k,Xvel)-velArr(i-1,j,k,Xvel))/(2.0*dx[0]);
                dvdy = (velArr(i, j+1, k,Yvel)-velArr(i,j-1,k,Yvel))/(2.0*dx[1]);
                divuArr(i,j,k,0) = dudx + dvdy;
            }
        }
    }
}

void
Projamrex::computeMagVorticity(const FArrayBox& statein, FArrayBox& der_out, const Box& bx, const Real* dx)
{
    Array4<Real const> const& vel = statein.array();
    Array4<Real> const& deriveArr = der_out.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    double vx, uy;
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                vx = (vel(i+1,j,k,Yvel)-vel(i-1,j,k,Yvel))/(2.0*dx[0]);
                uy = (vel(i,j+1,k,Xvel)-vel(i,j-1,k,Xvel))/(2.0*dx[1]);
                deriveArr(i,j,k,magvort) = (vx-uy);
            }
        }
    }
}

void
Projamrex::updateDensity(const FArrayBox& flux, FArrayBox& S_fab, const Box& bx, const Real* dx,
    const Real dt, int dim)
{
    Array4<Real const> const& fluxArray = flux.array();
    Array4<Real> const& stateNew = S_fab.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);

        for(int k = lo.z; k <= hi.z; k++)
        {
            for(int j = lo.y; j <= hi.y; j++)
            {
                for(int i = lo.x; i <= hi.x; i++)
                {
                    if(dim == 0)
                    {
                        stateNew(i,j,k,Density) = stateNew(i,j,k,Density) +
                            (dt/dx[0])*(fluxArray(i,j,k,Density)-fluxArray(i+1,j,k,Density));
                    }
                    else if(dim == 1)
                    {
                        stateNew(i,j,k,Density) = stateNew(i,j,k,Density) +
                            (dt/dx[1])*(fluxArray(i,j,k,0)-fluxArray(i,j+1,k,0));
                    }

                }
            }
        }
}

Real
Projamrex::upwind(const Real uL, const Real uR, const Real vel)
{
    if(vel > 0.0)
    {
        return uL;

    }
    else if(vel < 0.0)
    {
        return uR;
    }
    else
    {
        return 0.5*(uL+uR);
    }
}

void
Projamrex::calcTimeStep(const FArrayBox& statein, const Real* dx,
    amrex_real& dt)
{
    Array4<Real const> const& stateArray = statein.array();
    Dim3 lo = lbound(stateArray);
    Dim3 hi = ubound(stateArray);
    int Nv = stateArray.nComp();
    std::vector<Real> WCell(Nv);
    std::vector<Real> UCell(Nv);

    Real Cv = 0.0;
    Real Cc = 0.0;
    Real rho, u, v;
    Real advSpeed = 1.0;
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                u = stateArray(i,j,k,Xvel);
                v = stateArray(i,j,k,Yvel);
                rho = stateArray(i,j,k,Density);

                Cc = fmax(Cc, (fabs(u)/dx[0]+ fabs(v)/dx[1]));
                Cv = fmax(Cv, (mu/rho)*(1.0/(dx[0]*dx[0]) + 1.0/(dx[1]*dx[1])));
            }
        }
    }
    Cv = Cv*2.0;
     // dt = CFL*(dx[0]);
    dt = cfl/(Cv+Cc);
    if(dt <= 1.0e-15)
    {
        amrex::Abort("Time step too small!");
    }
    // std::cout << dt << std::endl;
}


//Estimate time step.
//
Real
Projamrex::estTimeStep (Real)
{
    // This is just a dummy value to start with
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    // std::cout << dx[1] << std::endl;
    const Real* prob_lo = geom.ProbLo();
    const Real cur_time = state[State_Type].curTime();
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab SBorder(grids, dmap, S_new.nComp(), NUM_GROW);
    MultiFab::Copy(SBorder, S_new, 0, 0, S_new.nComp(), S_new.nGrow());
    FillPatch(*this, SBorder, NUM_GROW, cur_time, State_Type, 0, S_new.nComp());
    Real dtCycle = 0.0;
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	FArrayBox uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
        const Box& bx = mfi.tilebox();

        const FArrayBox& statein = S_new[mfi];
        calcTimeStep(statein, dx, dtCycle);
        dt_est = fmin(dt_est, dtCycle);

	}
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    // dt_est *= cfl;

    if (verbose) {
	amrex::Print() << "Projamrex::estTimeStep at level " << level
                       << ":  dt_est = " << dt_est << std::endl;
    }

    return dt_est;
}

//
//Compute initial time step.
//
Real
Projamrex::initialTimeStep ()
{
    return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
Projamrex::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Compute new `dt'.
//
void
Projamrex::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        Projamrex& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1)
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    }
    else
    {
	//
	// Limit dt's by change_max * old dt
	//
	static Real change_max = 1.1;
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	}
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Do work after timestep().
//
void
Projamrex::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

#ifdef AMREX_PARTICLES
    if (TracerPC)
      {
        const int ncycle = parent->nCycle(level);

        if (iteration < ncycle || level == 0)
	  {
            int ngrow = (level == 0) ? 0 : iteration;

	    TracerPC->Redistribute(level, TracerPC->finestLevel(), ngrow);
	  }
      }
#endif
}

//
//Do work after regrid().
//
void
Projamrex::post_regrid (int lbase, int new_finest) {
#ifdef AMREX_PARTICLES
  if (TracerPC && level == lbase) {
      TracerPC->Redistribute(lbase);
  }
#endif
}

//
//Do work after a restart().
//
void
Projamrex::post_restart()
{
#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      BL_ASSERT(TracerPC == 0);
      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->Restart(parent->theRestartFile(), "Tracer");
    }
#endif
}

//
//Do work after init().
//
void
Projamrex::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//
void
Projamrex::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();
    Real rhoGradMaxCycle = 0.0;
    Real rhoGradMax = 0.0;
    const Real cur_time = state[State_Type].curTime();
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab SBorder(grids, dmap, S_new.nGrow(), NUM_GROW);
    MultiFab::Copy(SBorder, S_new, 0, 0, S_new.nGrow(), S_new.nGrow());
    FillPatch(*this, SBorder, NUM_GROW, cur_time, State_Type, 0, S_new.nGrow());
    MultiFab rhoGradFab(grids, dmap, 1, NUM_GROW);


#ifdef _OPENMP
#pragma omp parallel
// #pragma omp parallel reduction(min:dt_est)
#endif
    {
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
        const Box& bx = mfi.tilebox();
        const FArrayBox& statein = SBorder[mfi];
        FArrayBox& rhoGradIn = rhoGradFab[mfi];
        computeMaxGrad(statein, rhoGradIn, bx, dx, rhoGradMaxCycle);
        rhoGradMax = fmax(rhoGradMax, rhoGradMaxCycle);
	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
    }
    }
    ParallelDescriptor::ReduceRealMax(rhoGradMax);
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
                Vector<int>  itags;

        const FArrayBox& statein = SBorder[mfi];
        // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
        // So we are going to get a temporary integer array.
        const FArrayBox& rhoGradIn = rhoGradFab[mfi];

        const Box&  tilebx  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];
	    tagfab.get_itags(itags, tilebx);

        Array4<char> tagfabArray = tagfab.array();
            // data pointer and index space

	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebx.loVect();
	    const int*  thi     = tilebx.hiVect();

        performTaggingRhoGrad(tagfabArray,  tilebx, statein, rhoGradIn,rhoGradMax, dx);
	    // state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
		// 	BL_TO_FORTRAN_3D(S_new[mfi]),
		// 	&tagval, &clearval,
		// 	AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()),
		// 	AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);

	    // Now update the tags in the TagBox.
	    //
	    // tagfab.tags_and_untags(itags, tilebx);
	}
    }
}

void
Projamrex::computeMaxGrad(const FArrayBox& statein, FArrayBox& rhoGradIn, const Box& bx,
    const Real* dx, Real& rhoGradMax)
{
    Array4<Real const> const& stateArray = statein.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    Real rhoL, rhoC, rhoR, rhoU, rhoD;
    Real dRho_dx, dRho_dy;
    Real temp = 0.0;
    Array4<Real> const& rhoGrad = rhoGradIn.array();

    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {
                rhoL = stateArray(i,j+1,k,0);
                rhoC = stateArray(i+1, j+1, k, 0);
                rhoR = stateArray(i+2, j+1, k, 0);
                rhoU = stateArray(i+1, j+2, k, 0);
                rhoD = stateArray(i+1, j, k, 0);
                dRho_dx = (rhoR-rhoL)/(2*dx[0]);
                dRho_dy = (rhoU-rhoD)/(2*dx[1]);

                rhoGrad(i,j,k,0) = sqrt(dRho_dx*dRho_dx + dRho_dy*dRho_dy);

                if(temp < rhoGrad(i,j,k,0))
                {
                    temp = rhoGrad(i,j,k,0);
                }
            }
        }
    }
    // std::cout << "Break" << " ";

    rhoGradMax = temp;

}

void
Projamrex::performTaggingRhoGrad(Array4<char> tagfabArray, const Box& tilebx,
    const FArrayBox& statein,  const FArrayBox& rhoGradIn, Real rhoGradMax, const Real* dx)
{
    Array4<Real const> const& stateArray = statein.array();
    Dim3 lo = lbound(tilebx);
    Dim3 hi = ubound(tilebx);
    Real rhoL, rhoC, rhoR, rhoU, rhoD;
    Real dRho_dx, dRho_dy;
    Array4<Real const> const& rhoGrad = rhoGradIn.array();

    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {

                if(rhoGrad(i,j,k,0) >= 0.4*rhoGradMax)
                {
                    tagfabArray(i,j,k,0) = TagBox::SET;
                }
                else
                {
                    tagfabArray(i,j,k,0) = TagBox::CLEAR;

                }
            }
        }
    }

}

void
Projamrex::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);

    Geometry const* gg = AMReX::top()->getDefaultGeometry();

    // This tutorial code only supports Cartesian coordinates.
    if (! gg->IsCartesian()) {
	amrex::Abort("Please set geom.coord_sys = 0");
    }

    // This tutorial code only supports periodic boundaries.
    // if (! gg->isAllPeriodic()) {
	// amrex::Abort("Please set geom.is_periodic = 1 1 1");
    // }

#ifdef AMREX_PARTICLES
    pp.query("do_tracers", do_tracers);
#endif

    //
    // read tagging parameters from probin file
    //

    std::string probin_file("probin");

    ParmParse ppa("amr");
    ppa.query("probin_file",probin_file);

    int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

    // use a fortran routine to
    // read in tagging parameters from probin file
    // get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

}

void
Projamrex::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(State_Type),1.0,0,0,NUM_STATE,geom);

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = amrex::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        amrex::Print() << "Projamrex::reflux() at level " << level
                       << " : time = " << end << std::endl;
    }
}

void
Projamrex::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(State_Type);
}

void
Projamrex::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    Projamrex& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);

    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

#ifdef AMREX_PARTICLES
void
Projamrex::init_particles ()
{
  if (do_tracers and level == 0)
    {
      BL_ASSERT(TracerPC == nullptr);

      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->do_tiling = true;
      TracerPC->tile_size = IntVect(AMREX_D_DECL(1024000,4,4));

      AmrTracerParticleContainer::ParticleInitData pdata = {AMREX_D_DECL(0.0, 0.0, 0.0)};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);

      TracerPC->Redistribute();
    }
}
#endif

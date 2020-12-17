#ifndef _Projamrex_H_
#define _Projamrex_H_

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticles.H>
#endif

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <AMReX.H>
// #include <.H>
#include <memory>
#include <iostream>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_NodalProjector.H>
#ifdef _OPENMP
#include <omp.h>
#endif

enum StateType {State_Type = 0, Press_Type};

//
// AmrLevel-derived class for hyperbolic conservation equations
//

class Projamrex
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Projamrex ();
    //
    //The basic constructor.
    //
    Projamrex (amrex::Amr&     papa,
	         int             lev,
                 const amrex::Geometry& level_geom,
                 const amrex::BoxArray& bl,
                 const amrex::DistributionMapping& dm,
                 amrex::Real            time);
    //
    //The destructor.
    //
    virtual ~Projamrex () override;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (amrex::Amr&   papa,
                          std::istream& is,
			  bool          bReadSpecial = false) override;

    virtual void checkPoint (const std::string& dir,
			     std::ostream&      os,
			     amrex::VisMF::How  how = amrex::VisMF::NFiles,
			     bool               dump_old = true) override;

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData () override;
    //
    //Initialize data on this level from another Projamrex (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;

    void initBCs ();
    //
    void initialState(amrex::FArrayBox& statein, const amrex::Box& bx,
        amrex::RealBox gridloc, const amrex::Real* dx, int testNumber);

    void setFaceVelocity(amrex::FArrayBox& vel_x, amrex::FArrayBox& vel_y,
        const amrex::FArrayBox& statein, const amrex::Box& bx);

    void scalarAdvection(const amrex::FArrayBox& statein, amrex::FArrayBox& stateout, amrex::FArrayBox& flux_x,
    amrex::FArrayBox& flux_y, const amrex::FArrayBox& vel_x, const amrex::FArrayBox& vel_y, const amrex::Box& bx,
        const amrex::Real* dx, const amrex::Real dt);

    void computeNormalAdvectionVelocity(const amrex::FArrayBox& statein, const amrex::FArrayBox& vel_x,
        const amrex::FArrayBox& vel_y, amrex::FArrayBox& advVelx, amrex::FArrayBox& advVely,
            const amrex::Box& bx, const amrex::Real* dx, const amrex::Real dt);

    void computeTangentialAdvectionVelocity(const amrex::FArrayBox& statein, const amrex::FArrayBox& vel_x,
    const amrex::FArrayBox& vel_y, amrex::FArrayBox& advVelx, amrex::FArrayBox& advVely,
        const amrex::Box& bx, const amrex::Real* dx, const amrex::Real dt);

    void advectVelocity(const amrex::FArrayBox& advVelCC, const amrex::FArrayBox& u_halfx,
    const amrex::FArrayBox& v_halfy, const amrex::FArrayBox& u_halfy, const amrex::FArrayBox& v_halfx,
    const amrex::FArrayBox& statein, amrex::FArrayBox& stateout,
        const amrex::Box& bx, const amrex::Real* dx, const amrex::Real dt);

    void advectScalars(const amrex::FArrayBox& statein, amrex::FArrayBox& stateout, const amrex::FArrayBox& u_edge,
    const amrex::FArrayBox& v_edge, const amrex::FArrayBox& u_half, const amrex::FArrayBox& v_half, amrex::FArrayBox& flux_x,
    amrex::FArrayBox& flux_y, const amrex::Box& bx, const amrex::Real* dx, const amrex::Real dt);


    void velocityAdvection(const amrex::FArrayBox& statein, amrex::FArrayBox& stateout, amrex::FArrayBox& flux_x,
     amrex::FArrayBox& flux_y, const amrex::FArrayBox& vel_x, const amrex::FArrayBox& vel_y,
     const amrex::Box& bx, const amrex::Real* dx, const amrex::Real dt);

     void velocityDiffusion(const amrex::FArrayBox& statein, amrex::FArrayBox& stateout, const amrex::Box& bx,
        const amrex::Real* dx, const amrex::Real dt);

    void computeAdvectionTerm(amrex::FArrayBox& vel_x, amrex::FArrayBox& vel_y,
        amrex::FArrayBox& adv_x, amrex::FArrayBox& adv_y,
            const amrex::Box& bx, const amrex::Real* dx);

    void computeTemporaryVelocity(amrex::FArrayBox& adv_x, amrex::FArrayBox& adv_y,
        amrex::FArrayBox& tempvel_x, amrex::FArrayBox& tempvel_y,
            amrex::FArrayBox& vel_x, amrex::FArrayBox& vel_y,const amrex::Box& bx,
            const amrex::Real* dx, const amrex::Real dt);

    void computeAdvectiveFlux(const amrex::FArrayBox& vel, amrex::FArrayBox& rho,
                 amrex::FArrayBox& flux, const amrex::Box& bx,
                 const amrex::Real* dx, int dim);

    void computeDivergence(amrex::FArrayBox& divu, const amrex::FArrayBox& vel,
        const amrex::Box& bx, const amrex::Real* dx);

    void updateDensity(const amrex::FArrayBox& flux, amrex::FArrayBox& S_fab,
        const amrex::Box& bx, const amrex::Real* dx, const amrex::Real dt, int dim);

    void computeMagVorticity(const amrex::FArrayBox& statein, amrex::FArrayBox& der_out,
        const amrex::Box& bx, const amrex::Real* dx);


    void doNodalProjection(amrex::MultiFab& vel, amrex::MultiFab& pressure,
        amrex::MultiFab& S_cc, amrex::MultiFab& rho, amrex::Real dt);

    void doMacProjection(amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& vel);

    void fillphysbc_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    fetchBCsForProjection (amrex::Orientation::Side side) const noexcept;

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    fetchBCsForMacProjection (amrex::Orientation::Side side) const noexcept;

    amrex::Real upwind(const amrex::Real uL, const amrex::Real uR, const amrex::Real vel);

    void calcTimeStep(const amrex::FArrayBox& statein, const amrex::Real* dx,
        amrex_real& dt);

    void computeMaxGrad(const amrex::FArrayBox& statein, amrex::FArrayBox& rhoGradIn,
    const amrex::Box& bx, const amrex::Real* dx, amrex::Real& rhoGradMax);

    void performTaggingRhoGrad(amrex::Array4<char> tagfabArray, const amrex::Box& tilebx,
        const amrex::FArrayBox& statein,  const amrex::FArrayBox& rhoGradIn,
            amrex::Real rhoGradMax, const amrex::Real* dx);
    //Advance grids at this level in time.
    //
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;
    //
    //Estimate time step.
    //
    amrex::Real estTimeStep (amrex::Real dt_old);
    //
    //Compute initial time step.
    //
    amrex::Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) override;

    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase, int new_finest) override;
    //
    //Do work after a restart().
    //
    virtual void post_restart () override;
    //
    //Do work after init().
    //
    virtual void post_init (amrex::Real stop_time) override;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real         time,
			   int          n_error_buf = 0, int ngrow = 0) override;

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* theTracerPC () { return TracerPC.get(); }
#endif

    static int  NUM_STATE;
    static int  NUM_GROW;
    static amrex::BCRec       phys_bc;

protected:

    static void read_params ();

    // inline functioned implemented below
    Projamrex& getLevel (int lev);
    amrex::FluxRegister& getFluxReg ();
    amrex::FluxRegister& getFluxReg (int lev);

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

    //
    // The data.
    //
    amrex::FluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static int          verbose;
    static amrex::Real  cfl;
    static amrex::Real  mu;
    static int          do_reflux;
    static int          num_state_type;
    static int          m_ntrac;
    static int          testNumber;

    enum struct BC
    {
        noSlipWall, slipWall, periodic, undefined
    };

     amrex::Array<BC, 2*AMREX_SPACEDIM> m_bc_type;
     amrex::Array<amrex::Real, 2*AMREX_SPACEDIM> m_bc_press;
     amrex::Array<amrex::Real, 2*AMREX_SPACEDIM> m_bc_den;
     amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM>,
                            2*AMREX_SPACEDIM> m_bc_vel;
     amrex::Array<amrex::Vector<amrex::Real> , AMREX_SPACEDIM*2> m_bc_trac;

     amrex::Vector<amrex::BCRec> m_bcrec_vel;
     amrex::Vector<amrex::BCRec> m_bcrec_den;
     amrex::Vector<amrex::BCRec> m_bcrec_trac;
     amrex::Vector<amrex::BCRec> m_bcrec_press;

#ifdef AMREX_PARTICLES
    void init_particles ();
    static int       do_tracers;
    static std::unique_ptr<amrex::AmrTracerParticleContainer> TracerPC;
#endif
};

//
// Inlines.
//

inline
Projamrex&
Projamrex::getLevel (int lev)
{
    return *(Projamrex *) &parent->getLevel(lev);
}

inline
amrex::FluxRegister&
Projamrex::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
amrex::FluxRegister&
Projamrex::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}
#endif /*_Projamrex_H_*/

#ifndef BCFUNCS_H
#define BCFUNCS_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

struct ProjVelFill
{

  int probtype;
  amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM*2> bcv;

  AMREX_GPU_HOST
  constexpr ProjVelFill (int a_probtype,
                           amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM>,
                                           AMREX_SPACEDIM*2> const& a_bcv)
      : probtype(a_probtype), bcv(a_bcv) {}

  AMREX_GPU_DEVICE
  void operator() (const IntVect& iv, Array4<Real> const& vel,
                   const int dcomp, const int numcomp,
                   GeometryData const& geom, const Real time,
                   const BCRec* bcr, const int bcomp,
                   const int orig_comp) const
  {


    using namespace amrex;

    // do something for external Dirichlet (BCType::ext_dir)
      const int i = iv[0];
      const int j = iv[1];
      const int k = iv[2];

      const Box& dombox = geom.Domain();

      // external Dirichlet for cell iv
      for(int vcomp = 0; vcomp < AMREX_SPACEDIM; ++vcomp)
      {
          const BCRec& bc = bcr[bcomp+vcomp];
          if(bc.lo(0) == BCType::ext_dir and i < dombox.smallEnd(0))
          {
            vel(i,j,k,vcomp) = bcv[Orientation(Direction::x,Orientation::low)][vcomp];
          }
          else if(bc.hi(0) == BCType::ext_dir and i > dombox.bigEnd(0));
          {
            vel(i,j,k,vcomp) = bcv[Orientation(Direction::x,Orientation::high)][vcomp];
          }
          if(bc.lo(1) == BCType::ext_dir and j < dombox.smallEnd(1))
          {
            vel(i,j,k,vcomp) = bcv[Orientation(Direction::y,Orientation::low)][vcomp];
          }
          else if(bc.hi(1) == BCType::ext_dir and j > dombox.bigEnd(1));
          {
            vel(i,j,k,vcomp) = bcv[Orientation(Direction::y,Orientation::high)][vcomp];
          }

      }
  }
};

struct ProjDenFill
{
    int probtype;
    amrex::Array<amrex::Real,AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr ProjDenFill (int a_probtype, amrex::Array<amrex::Real, AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& rho,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const Box& domain_box = geom.Domain();
        const BCRec& bc = bcr[bcomp];

        if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
        {
            rho(i,j,k) = bcv[Orientation(Direction::x,Orientation::low)];
        }
        else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
        {
            rho(i,j,k) = bcv[Orientation(Direction::x,Orientation::high)];
        }

        if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
        {
            rho(i,j,k) = bcv[Orientation(Direction::y,Orientation::low)];
        }
        else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
        {
            rho(i,j,k) = bcv[Orientation(Direction::y,Orientation::high)];
        }

    }
};

struct ProjPressFill
{
    int probtype;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr ProjPressFill (int a_probtype, amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& press,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const Box& domain_box = geom.Domain();
        const BCRec& bc = bcr[bcomp];

        if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
        {
            press(i,j,k) = bcv[Orientation(Direction::x,Orientation::low)];
        }
        else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
        {
            press(i,j,k) = bcv[Orientation(Direction::x,Orientation::high)];
        }

        if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
        {
            press(i,j,k) = bcv[Orientation(Direction::y,Orientation::low)];
        }
        else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
        {
            press(i,j,k) = bcv[Orientation(Direction::y,Orientation::high)];
        }

    }
};

struct ProjTracFill
{
    int probtype, ntrac;
    amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr ProjTracFill (int a_probtype, int a_ntrac,
                              amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), ntrac(a_ntrac), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& tracer,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const Box& domain_box = geom.Domain();

        for (int n = 0; n < ntrac; ++n)
        {
            const BCRec& bc = bcr[bcomp+n];

            if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::x,Orientation::low)][n];
            }
            else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::x,Orientation::high)][n];
            }

            if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::y,Orientation::low)][n];
            }
            else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::y,Orientation::high)][n];
            }

        }
    }
};

#endif
